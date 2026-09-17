#!/usr/bin/env python3
"""Build leakage-clean EVN/SVN chromosome datasets from genomic backgrounds.

Exact genomic-background loci are deduplicated globally in split-priority order:
chr1 validation, chr2 test, then chr3--chr22 and chrX training. SVN retains all
eligible unique genomic-background loci and defines the common per-class target. EVN is then
allocated to that target with Hamilton chromosome-proportional quotas.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import tempfile
from pathlib import Path
from typing import TextIO

from grch38_training_utils import (
    SourcePair,
    copy_selected_pair,
    count_bed,
    count_fasta,
    sample_candidates,
    unique_coordinate_candidates,
    validate_unique_bed,
    validate_unique_fasta,
)


CHROMS = tuple(str(i) for i in range(1, 23)) + ("X",)
TRAIN_CHROMS = tuple(c for c in CHROMS if c not in ("1", "2"))


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def split_for(chrom: str) -> str:
    return "validation" if chrom == "1" else "test" if chrom == "2" else "train"


def source(repo: Path, dataset: str, chrom: str, cls: str) -> SourcePair:
    directory = repo / f"data/processed/{dataset}/GRCh38/chroms/chrom{chrom}"
    suffix = "ELS" if dataset == "evn" else "Silencer_ELS"
    stem = f"chrom{chrom}_GRCh38_{suffix}"
    prefix = "" if cls == "pos" else "neg1x_"
    bed = directory / f"{prefix}{stem}.bed"
    fasta = directory / f"{prefix}{stem}.fa"
    n_bed, n_fasta = count_bed(bed), count_fasta(fasta)
    if n_bed != n_fasta:
        raise ValueError(f"{dataset} chr{chrom} {cls}: BED/FASTA mismatch")
    return SourcePair(bed, fasta, n_bed)


def hamilton(counts: dict[str, int], target: int) -> dict[str, int]:
    total = sum(counts.values())
    exact = {c: target * counts[c] / total for c in CHROMS}
    quotas = {c: int(exact[c]) for c in CHROMS}
    remainder = target - sum(quotas.values())
    order = {c: i for i, c in enumerate(CHROMS)}
    ranked = sorted(CHROMS, key=lambda c: (-(exact[c] - quotas[c]), order[c]))
    for chrom in ranked[:remainder]:
        quotas[chrom] += 1
    if sum(quotas.values()) != target:
        raise AssertionError("Hamilton quotas do not sum to target")
    return quotas


def plan(repo: Path, seed: int) -> tuple[dict, int]:
    planned: dict[str, dict[str, dict]] = {"svn": {}, "evn": {}}
    svn_seen: set[tuple[str, int, int]] = set()
    svn_total = 0
    for ci, chrom in enumerate(CHROMS):
        pos, neg = source(repo, "svn", chrom, "pos"), source(repo, "svn", chrom, "neg")
        candidates = unique_coordinate_candidates(neg, svn_seen)
        quota = min(pos.count, len(candidates))
        pos_indices = set(random.Random(seed + ci).sample(range(pos.count), quota))
        neg_indices, coords = sample_candidates(candidates, quota, seed + 100_000 + ci)
        svn_seen.update(coords)
        planned["svn"][chrom] = {
            "pos": pos, "neg": neg, "pos_indices": pos_indices,
            "neg_indices": neg_indices, "quota": quota,
            "unique_negative_eligible": len(candidates),
        }
        svn_total += quota

    evn_sources = {
        chrom: {
            "pos": source(repo, "evn", chrom, "pos"),
            "neg": source(repo, "evn", chrom, "neg"),
        }
        for chrom in CHROMS
    }
    quotas = hamilton({c: evn_sources[c]["pos"].count for c in CHROMS}, svn_total)
    evn_seen: set[tuple[str, int, int]] = set()
    for ci, chrom in enumerate(CHROMS):
        pos, neg = evn_sources[chrom]["pos"], evn_sources[chrom]["neg"]
        candidates = unique_coordinate_candidates(neg, evn_seen)
        quota = quotas[chrom]
        if quota > len(candidates):
            raise ValueError(
                f"EVN chr{chrom}: quota {quota:,} exceeds {len(candidates):,} unique real negatives"
            )
        pos_indices = set(random.Random(seed + 200_000 + ci).sample(range(pos.count), quota))
        neg_indices, coords = sample_candidates(candidates, quota, seed + 300_000 + ci)
        evn_seen.update(coords)
        planned["evn"][chrom] = {
            "pos": pos, "neg": neg, "pos_indices": pos_indices,
            "neg_indices": neg_indices, "quota": quota,
            "unique_negative_eligible": len(candidates),
        }
    return planned, svn_total


def output_paths(root: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {
        "manifest": root / "chromosome_training_manifest.tsv",
        "summary": root / "chromosome_training_summary.json",
    }
    for split in ("validation", "test", "train_100"):
        for cls in ("pos", "neg"):
            for ext in ("bed", "fa"):
                paths[f"{split}_{cls}_{ext}"] = root / f"{split}_{cls}.{ext}"
    for chrom in TRAIN_CHROMS:
        for cls in ("pos", "neg"):
            for ext in ("bed", "fa"):
                paths[f"chr{chrom}_{cls}_{ext}"] = (
                    root / "chromosomes" / f"chrom{chrom}" / f"train_{cls}.{ext}"
                )
    return paths


def open_outputs(paths: dict[str, Path], overwrite: bool) -> tuple[dict[str, Path], dict[str, TextIO]]:
    existing = [p for p in paths.values() if p.exists()]
    if existing and not overwrite:
        raise FileExistsError("Outputs exist; pass --overwrite")
    temporary, handles = {}, {}
    for name, path in paths.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        handle = tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", newline="", dir=path.parent,
            prefix=f".{path.name}.", suffix=".tmp", delete=False,
        )
        temporary[name], handles[name] = Path(handle.name), handle
    return temporary, handles


def write_dataset(repo: Path, dataset: str, plan_data: dict[str, dict], seed: int, target: int, overwrite: bool) -> None:
    root = repo / f"data/processed/{dataset}/GRCh38/training"
    paths = output_paths(root)
    temporary, handles = open_outputs(paths, overwrite)
    manifest = csv.writer(handles["manifest"], delimiter="\t", lineterminator="\n")
    manifest.writerow((
        "chromosome", "split", "positive_available", "real_negative_available",
        "unique_real_negative_eligible", "selected_per_class",
        "positive_bed", "positive_fasta", "negative_bed", "negative_fasta",
    ))
    split_counts = {"train": 0, "validation": 0, "test": 0}
    try:
        for chrom in CHROMS:
            item = plan_data[chrom]
            split = split_for(chrom)
            manifest.writerow((
                f"chr{chrom}", split, item["pos"].count, item["neg"].count,
                item["unique_negative_eligible"], item["quota"], item["pos"].bed,
                item["pos"].fasta, item["neg"].bed, item["neg"].fasta,
            ))
            if split == "train":
                targets = [f"chr{chrom}", "train_100"]
            else:
                targets = [split]
            for cls in ("pos", "neg"):
                src = item[cls]
                indices = item[f"{cls}_indices"]
                for destination in targets:
                    copy_selected_pair(
                        src, indices,
                        handles[f"{destination}_{cls}_bed"],
                        handles[f"{destination}_{cls}_fa"],
                    )
            split_counts[split] += item["quota"]
        summary = {
            "dataset": dataset,
            "negative_kind": "genomic",
            "seed": seed,
            "common_target_per_class": target,
            "strategy": "global exact-locus deduplication; chr1 validation, chr2 test",
            "records_per_class": split_counts,
        }
        handles["summary"].write(json.dumps(summary, indent=2, sort_keys=True) + "\n")
        for handle in handles.values():
            handle.flush(); os.fsync(handle.fileno()); handle.close()
        for name, tmp in temporary.items():
            os.replace(tmp, paths[name])
    except BaseException:
        for handle in handles.values():
            if not handle.closed:
                handle.close()
        for tmp in temporary.values():
            tmp.unlink(missing_ok=True)
        raise

    for split, expected in (("validation", split_counts["validation"]), ("test", split_counts["test"]), ("train_100", split_counts["train"])):
        for cls in ("pos", "neg"):
            n_bed = count_bed(paths[f"{split}_{cls}_bed"])
            n_fasta = validate_unique_fasta(paths[f"{split}_{cls}_fa"])
            if n_bed != expected or n_fasta != expected:
                raise AssertionError(f"{dataset} {split} {cls}: output count mismatch")
    for cls in ("pos", "neg"):
        observed = validate_unique_bed(
            [paths[f"{split}_{cls}_bed"] for split in ("validation", "test", "train_100")],
            f"{dataset} genomic {cls}",
        )
        if observed != target:
            raise AssertionError(f"{dataset} {cls}: expected {target:,} unique loci, got {observed:,}")
    print(f"{dataset}: {target:,} records per class; validation={split_counts['validation']:,}, test={split_counts['test']:,}, train={split_counts['train']:,}")


def main() -> None:
    args = parse_args()
    repo = args.repo.resolve()
    planned, target = plan(repo, args.seed)
    for dataset in ("svn", "evn"):
        write_dataset(repo, dataset, planned[dataset], args.seed, target, args.overwrite)
    print("Leakage-clean real chromosome benchmark: PASS")


if __name__ == "__main__":
    main()
