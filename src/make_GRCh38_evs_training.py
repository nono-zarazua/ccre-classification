#!/usr/bin/env python3
"""Build a chromosome-matched enhancer-versus-silencer GRCh38 benchmark.

All cleaned SVN loci on chr1--chr22/X are retained. For each chromosome, the
clean EVN source is downsampled to the corresponding SVN count with a fixed
seed. Enhancers are the positive class and silencers are the negative class.
Chr1 is validation, chr2 is held-out test, and chr3--chr22/X are training.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import shutil
import subprocess
import tempfile
from pathlib import Path


CHROMOSOMES = tuple(str(i) for i in range(1, 23)) + ("X",)
TRAIN_CHROMOSOMES = tuple(c for c in CHROMOSOMES if c not in ("1", "2"))


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--reference", type=Path, default=Path("data/reference/hg38.fa"))
    parser.add_argument(
        "--bedtools",
        type=Path,
        default=Path("/home/zarazuanav/.conda/envs/bstool/bin/bedtools"),
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def source_bed(repo: Path, task: str, chromosome: str) -> Path:
    suffix = "ELS" if task == "evn" else "Silencer_ELS"
    return (
        repo
        / f"data/processed/{task}/GRCh38/chroms/chrom{chromosome}"
        / f"chrom{chromosome}_GRCh38_{suffix}.bed"
    )


def read_bed(path: Path, chromosome: str) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty BED: {path}")
    expected = f"chr{chromosome}"
    records: list[str] = []
    seen_ids: set[str] = set()
    seen_loci: set[tuple[str, int, int]] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"{path}:{line_number}: expected at least four columns")
            if fields[0] != expected:
                raise ValueError(
                    f"{path}:{line_number}: expected {expected}, observed {fields[0]}"
                )
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")
            if fields[3] in seen_ids:
                raise ValueError(f"{path}:{line_number}: duplicate ID {fields[3]}")
            locus = (fields[0], start, end)
            if locus in seen_loci:
                raise ValueError(f"{path}:{line_number}: duplicate locus {locus}")
            seen_ids.add(fields[3])
            seen_loci.add(locus)
            records.append("\t".join(fields) + "\n")
    return records


def assert_no_cross_class_overlap(
    enhancers: list[str], silencers: list[str], chromosome: str
) -> None:
    intervals = sorted(
        [(int(r.split("\t")[1]), int(r.split("\t")[2]), "enhancer") for r in enhancers]
        + [(int(r.split("\t")[1]), int(r.split("\t")[2]), "silencer") for r in silencers]
    )
    active: list[tuple[int, str]] = []
    for start, end, label in intervals:
        active = [(other_end, other_label) for other_end, other_label in active if other_end > start]
        if any(other_label != label for _, other_label in active):
            raise ValueError(f"chr{chromosome}: enhancer/silencer interval overlap remains")
        active.append((end, label))


def temporary_path(destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    )
    path = Path(handle.name)
    handle.close()
    return path


def write_lines_atomic(destination: Path, records: list[str]) -> None:
    temporary = temporary_path(destination)
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.writelines(records)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def extract_fasta(bedtools: Path, reference: Path, bed: Path, destination: Path) -> None:
    temporary = temporary_path(destination)
    try:
        subprocess.run(
            [
                str(bedtools), "getfasta", "-fi", str(reference), "-bed", str(bed),
                "-fo", str(temporary), "-name",
            ],
            check=True,
        )
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def count_fasta(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        return sum(line.startswith(">") for line in handle)


def fasta_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    current: int | None = None
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            text = line.strip()
            if text.startswith(">"):
                if current is not None:
                    lengths.append(current)
                current = 0
            elif text:
                if current is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                current += len(text)
    if current is not None:
        lengths.append(current)
    return lengths


def bed_lengths(path: Path) -> list[int]:
    lengths: list[int] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            lengths.append(int(fields[2]) - int(fields[1]))
    return lengths


def copy_atomic(source: Path, destination: Path) -> None:
    temporary = temporary_path(destination)
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def concatenate_atomic(sources: list[Path], destination: Path) -> None:
    temporary = temporary_path(destination)
    try:
        with temporary.open("wb") as output:
            for source in sources:
                with source.open("rb") as input_handle:
                    shutil.copyfileobj(input_handle, output)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> None:
    args = parse_args()
    repo = args.repo.resolve()
    reference = args.reference if args.reference.is_absolute() else repo / args.reference
    bedtools = args.bedtools.expanduser().resolve()
    if not reference.is_file() or reference.stat().st_size == 0:
        raise FileNotFoundError(f"Missing reference FASTA: {reference}")
    if not bedtools.is_file() or not os.access(bedtools, os.X_OK):
        raise FileNotFoundError(f"Missing executable bedtools: {bedtools}")

    root = repo / "data/processed/evs/GRCh38"
    chrom_root = root / "chroms"
    training_root = root / "real"
    selected: dict[str, dict[str, Path | int]] = {}
    manifest_rows: list[tuple[object, ...]] = []

    for index, chromosome in enumerate(CHROMOSOMES):
        evn_source = source_bed(repo, "evn", chromosome)
        svn_source = source_bed(repo, "svn", chromosome)
        enhancers_available = read_bed(evn_source, chromosome)
        silencers = read_bed(svn_source, chromosome)
        quota = len(silencers)
        if quota > len(enhancers_available):
            raise ValueError(
                f"chr{chromosome}: {quota:,} silencers exceed "
                f"{len(enhancers_available):,} available enhancers"
            )
        chosen = sorted(
            random.Random(args.seed + index).sample(range(len(enhancers_available)), quota)
        )
        enhancers = [enhancers_available[i] for i in chosen]
        assert_no_cross_class_overlap(enhancers, silencers, chromosome)

        directory = chrom_root / f"chrom{chromosome}"
        enhancer_bed = directory / f"chrom{chromosome}_GRCh38_enhancer.bed"
        enhancer_fasta = directory / f"chrom{chromosome}_GRCh38_enhancer.fasta"
        silencer_bed = directory / f"chrom{chromosome}_GRCh38_silencer.bed"
        silencer_fasta = directory / f"chrom{chromosome}_GRCh38_silencer.fasta"
        outputs = (enhancer_bed, enhancer_fasta, silencer_bed, silencer_fasta)
        if not args.overwrite and any(path.exists() for path in outputs):
            raise FileExistsError(f"chr{chromosome}: outputs exist; pass --overwrite")
        write_lines_atomic(enhancer_bed, enhancers)
        write_lines_atomic(silencer_bed, silencers)
        extract_fasta(bedtools, reference, enhancer_bed, enhancer_fasta)
        extract_fasta(bedtools, reference, silencer_bed, silencer_fasta)
        if count_fasta(enhancer_fasta) != quota or count_fasta(silencer_fasta) != quota:
            raise AssertionError(f"chr{chromosome}: BED/FASTA count mismatch")
        if fasta_lengths(enhancer_fasta) != bed_lengths(enhancer_bed):
            raise AssertionError(f"chr{chromosome}: enhancer BED/FASTA length mismatch")
        if fasta_lengths(silencer_fasta) != bed_lengths(silencer_bed):
            raise AssertionError(f"chr{chromosome}: silencer BED/FASTA length mismatch")

        selected[chromosome] = {
            "enhancer_bed": enhancer_bed,
            "enhancer_fasta": enhancer_fasta,
            "silencer_bed": silencer_bed,
            "silencer_fasta": silencer_fasta,
            "count": quota,
        }
        split = "validation" if chromosome == "1" else "test" if chromosome == "2" else "train"
        manifest_rows.append(
            (
                f"chr{chromosome}", split, len(enhancers_available), quota, quota,
                evn_source.relative_to(repo), svn_source.relative_to(repo),
                enhancer_bed.relative_to(repo), silencer_bed.relative_to(repo),
            )
        )

    for chromosome in TRAIN_CHROMOSOMES:
        shard = training_root / "chromosomes" / f"chrom{chromosome}"
        item = selected[chromosome]
        copy_atomic(item["enhancer_bed"], shard / "train_pos.bed")
        copy_atomic(item["enhancer_fasta"], shard / "train_pos.fa")
        copy_atomic(item["silencer_bed"], shard / "train_neg.bed")
        copy_atomic(item["silencer_fasta"], shard / "train_neg.fa")

    for chromosome, split in (("1", "validation"), ("2", "test")):
        item = selected[chromosome]
        copy_atomic(item["enhancer_bed"], training_root / f"{split}_pos.bed")
        copy_atomic(item["enhancer_fasta"], training_root / f"{split}_pos.fa")
        copy_atomic(item["silencer_bed"], training_root / f"{split}_neg.bed")
        copy_atomic(item["silencer_fasta"], training_root / f"{split}_neg.fa")

    concatenate_atomic(
        [selected[c]["enhancer_bed"] for c in TRAIN_CHROMOSOMES],
        training_root / "train_100_pos.bed",
    )
    concatenate_atomic(
        [selected[c]["enhancer_fasta"] for c in TRAIN_CHROMOSOMES],
        training_root / "train_100_pos.fa",
    )
    concatenate_atomic(
        [selected[c]["silencer_bed"] for c in TRAIN_CHROMOSOMES],
        training_root / "train_100_neg.bed",
    )
    concatenate_atomic(
        [selected[c]["silencer_fasta"] for c in TRAIN_CHROMOSOMES],
        training_root / "train_100_neg.fa",
    )

    manifest_path = root / "evs_manifest.tsv"
    temporary = temporary_path(manifest_path)
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "chromosome", "split", "enhancers_available", "enhancers_selected",
                "silencers_selected", "enhancer_source", "silencer_source",
                "enhancer_output", "silencer_output",
            )
        )
        writer.writerows(manifest_rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, manifest_path)

    counts_path = chrom_root / "chromosome_counts.csv"
    temporary = temporary_path(counts_path)
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            ("chromosome", "split", "enhancers_available", "enhancers_selected", "silencers")
        )
        for row in manifest_rows:
            writer.writerow(row[:5])
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, counts_path)

    split_counts = {
        "validation": int(selected["1"]["count"]),
        "test": int(selected["2"]["count"]),
        "training": sum(int(selected[c]["count"]) for c in TRAIN_CHROMOSOMES),
    }
    summary_path = root / "evs_summary.json"
    temporary = temporary_path(summary_path)
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "task": "enhancer-versus-silencer",
                "positive_class": "enhancer",
                "negative_class": "silencer",
                "seed": args.seed,
                "sampling": "all silencers; enhancers sampled to exact silencer count per chromosome",
                "excluded_chromosomes": ["chrY"],
                "split_chromosomes": {
                    "validation": ["chr1"],
                    "test": ["chr2"],
                    "training": [f"chr{c}" for c in TRAIN_CHROMOSOMES],
                },
                "records_per_class": split_counts,
                "total_per_class": sum(split_counts.values()),
            },
            handle,
            indent=2,
            sort_keys=True,
        )
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, summary_path)

    print(f"EVS validation per class: {split_counts['validation']:,}")
    print(f"EVS test per class:       {split_counts['test']:,}")
    print(f"EVS training per class:   {split_counts['training']:,}")
    print(f"EVS total per class:      {sum(split_counts.values()):,}")
    print("EVS chromosome matching and zero cross-class interval overlap: PASS")


if __name__ == "__main__":
    main()
