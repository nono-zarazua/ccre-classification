#!/usr/bin/env python3
"""Build a chromosome-proportional benchmark with fixed train/validation/test splits.

BED and FASTA files created by the existing chromosome pipeline correspond by
record position, not by a shared identifier. This script therefore samples a
record index once and applies it to both representations. chrY is excluded.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


CHROMOSOMES = tuple(str(value) for value in range(1, 23)) + ("X",)
FIXED_SPLIT_CHROMOSOME = {"validation": "1", "test": "2"}
TRAINING_CHROMOSOMES = tuple(
    chromosome for chromosome in CHROMOSOMES
    if chromosome not in FIXED_SPLIT_CHROMOSOME.values()
)
FOLD_CHROMOSOMES = {
    1: ("1", "7", "16", "X"),
    2: ("2", "12", "17", "19"),
    3: ("3", "8", "9", "18", "22"),
    4: ("6", "4", "11", "13", "21"),
    5: ("5", "10", "15", "14", "20"),
}
OUTER_FOLDS = {1: (2, 3, 4), 2: (1, 3, 4), 3: (1, 2, 4), 4: (1, 2, 3)}
VALIDATION_CHROMOSOME = {1: "6", 2: "6", 3: "6", 4: "3"}
CHROMOSOME_FOLD = {
    chromosome: fold
    for fold, chromosomes in FOLD_CHROMOSOMES.items()
    for chromosome in chromosomes
}


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    text: str


@dataclass(frozen=True)
class SourcePair:
    chromosome: str
    class_name: str
    bed_path: Path
    fasta_path: Path
    bed_records: tuple[str, ...]
    fasta_records: tuple[FastaRecord, ...]


def parse_args() -> argparse.Namespace:
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "Sample a fixed chromosome-proportional GRCh38 ELS benchmark, then "
            "route chr1 to validation, chr2 to test, and all others to training."
        )
    )
    parser.add_argument(
        "--chroms-dir",
        type=Path,
        default=repository / "data/processed/svn/GRCh38/chroms",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository / "data/processed/svn/GRCh38/training",
    )
    parser.add_argument(
        "--folds-dir",
        type=Path,
        default=repository / "data/processed/svn/GRCh38/folds",
    )
    parser.add_argument("--target-count", type=int, default=202_493)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing master outputs. Without this flag, fail safely.",
    )
    return parser.parse_args()


def source_paths(chroms_dir: Path, chromosome: str, class_name: str) -> tuple[Path, Path]:
    directory = chroms_dir / f"chrom{chromosome}"
    if class_name == "pos":
        stem = f"chrom{chromosome}_GRCh38_Silencer_ELS"
    else:
        stem = f"neg1x_chrom{chromosome}_GRCh38_Silencer_ELS"
    return directory / f"{stem}.bed", directory / f"{stem}.fa"


def read_bed(path: Path) -> tuple[str, ...]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"BED input is missing or empty: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        records = tuple(line if line.endswith("\n") else line + "\n" for line in handle)
    if any(not record.strip() for record in records):
        raise ValueError(f"Blank BED record found in {path}")
    return records


def read_fasta(path: Path) -> tuple[FastaRecord, ...]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"FASTA input is missing or empty: {path}")
    records: list[FastaRecord] = []
    identifier: str | None = None
    lines: list[str] = []
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if identifier is not None:
                    if len(lines) == 1:
                        raise ValueError(f"{path}: FASTA record {identifier!r} has no sequence")
                    records.append(FastaRecord(identifier, "".join(lines)))
                identifier = line[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                if identifier in seen:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA ID {identifier!r}")
                seen.add(identifier)
                lines = [line]
            else:
                if identifier is None:
                    if line.strip():
                        raise ValueError(f"{path}:{line_number}: sequence before first header")
                    continue
                lines.append(line)
    if identifier is not None:
        if len(lines) == 1:
            raise ValueError(f"{path}: FASTA record {identifier!r} has no sequence")
        text = "".join(lines)
        if not text.endswith(("\n", "\r")):
            text += "\n"
        records.append(FastaRecord(identifier, text))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return tuple(records)


def load_source(chroms_dir: Path, chromosome: str, class_name: str) -> SourcePair:
    bed_path, fasta_path = source_paths(chroms_dir, chromosome, class_name)
    bed_records = read_bed(bed_path)
    fasta_records = read_fasta(fasta_path)
    if len(bed_records) != len(fasta_records):
        raise ValueError(
            f"Positional BED/FASTA mismatch for chr{chromosome} {class_name}: "
            f"{len(bed_records):,} BED versus {len(fasta_records):,} FASTA records"
        )
    return SourcePair(
        chromosome, class_name, bed_path, fasta_path, bed_records, fasta_records
    )


def hamilton_quotas(counts: dict[str, int], target: int) -> dict[str, int]:
    """Allocate an exact target using proportional largest-remainder rounding."""
    total = sum(counts.values())
    if target <= 0:
        raise ValueError("--target-count must be positive")
    if target > total:
        raise ValueError(f"Target {target:,} exceeds {total:,} available positives")
    exact = {chrom: target * count / total for chrom, count in counts.items()}
    quotas = {chrom: int(value) for chrom, value in exact.items()}
    remaining = target - sum(quotas.values())
    chromosome_order = {chrom: index for index, chrom in enumerate(CHROMOSOMES)}
    ranked = sorted(
        counts,
        key=lambda chrom: (-(exact[chrom] - quotas[chrom]), chromosome_order[chrom]),
    )
    for chromosome in ranked[:remaining]:
        quotas[chromosome] += 1
    if sum(quotas.values()) != target:
        raise AssertionError("Hamilton allocation did not sum to the requested target")
    return quotas


def sample_indices(size: int, quota: int, seed: int) -> tuple[int, ...]:
    if quota > size:
        raise ValueError(f"Requested {quota:,} records from a pool of {size:,}")
    # Sorting preserves source order in output while membership remains randomized.
    return tuple(sorted(random.Random(seed).sample(range(size), quota)))


def output_paths(output_dir: Path, folds_dir: Path) -> dict[str, Path]:
    paths = {
        "pos_bed": output_dir / "train_100_pos.bed",
        "pos_fa": output_dir / "train_100_pos.fa",
        "neg_bed": output_dir / "train_100_neg.bed",
        "neg_fa": output_dir / "train_100_neg.fa",
        "manifest": output_dir / "train_100_manifest.tsv",
        "summary": output_dir / "train_100_summary.json",
        "validation_pos_bed": output_dir / "validation_pos.bed",
        "validation_pos_fa": output_dir / "validation_pos.fa",
        "validation_neg_bed": output_dir / "validation_neg.bed",
        "validation_neg_fa": output_dir / "validation_neg.fa",
        "test_pos_bed": output_dir / "test_pos.bed",
        "test_pos_fa": output_dir / "test_pos.fa",
        "test_neg_bed": output_dir / "test_neg.bed",
        "test_neg_fa": output_dir / "test_neg.fa",
    }
    for fold in FOLD_CHROMOSOMES:
        fold_dir = folds_dir / f"fold{fold}"
        paths[f"fold{fold}_pos_bed"] = fold_dir / f"fold{fold}.bed"
        paths[f"fold{fold}_pos_fa"] = fold_dir / f"fold{fold}.fa"
        paths[f"fold{fold}_neg_bed"] = fold_dir / f"neg1x_fold{fold}.bed"
        paths[f"fold{fold}_neg_fa"] = fold_dir / f"neg1x_fold{fold}.fa"
    for outer in OUTER_FOLDS:
        outer_dir = output_dir / f"outer{outer}"
        for split in ("train", "validation"):
            for class_name in ("pos", "neg"):
                for extension in ("bed", "fa"):
                    paths[f"outer{outer}_{split}_{class_name}_{extension}"] = (
                        outer_dir / f"{split}_{class_name}.{extension}"
                    )
    return paths


def record_destinations(chromosome: str) -> tuple[str, ...]:
    """Return fold and outer destinations for a sampled source-chromosome record."""
    fold = CHROMOSOME_FOLD[chromosome]
    destinations = [f"fold{fold}"]
    for outer, included_folds in OUTER_FOLDS.items():
        if chromosome == VALIDATION_CHROMOSOME[outer]:
            destinations.append(f"outer{outer}_validation")
        elif fold in included_folds:
            destinations.append(f"outer{outer}_train")
    return tuple(destinations)


def bed_coordinates(record: str, path: Path, source_index: int) -> tuple[str, str, str]:
    fields = record.rstrip("\r\n").split("\t")
    if len(fields) < 3:
        raise ValueError(f"{path}: record {source_index + 1} has fewer than 3 columns")
    return fields[0], fields[1], fields[2]


def create_temporary(path: Path) -> tuple[Path, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    return Path(handle.name), handle


def write_outputs(
    paths: dict[str, Path],
    sources: dict[str, dict[str, SourcePair]],
    selections: dict[str, dict[str, tuple[int, ...]]],
    quotas: dict[str, int],
    target: int,
    seed: int,
    overwrite: bool,
) -> None:
    existing = [path for path in paths.values() if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing master output(s):\n  "
            + "\n  ".join(str(path) for path in existing)
        )

    temporary_paths: dict[str, Path] = {}
    handles: dict[str, object] = {}
    try:
        for name, path in paths.items():
            temporary_paths[name], handles[name] = create_temporary(path)

        manifest = csv.writer(handles["manifest"], delimiter="\t", lineterminator="\n")
        manifest.writerow(
            ("class", "fixed_split", "source_chromosome", "assigned_fold",
             "source_record_1based", "bed_chromosome", "bed_start", "bed_end",
             "fasta_id")
        )

        for class_name in ("pos", "neg"):
            for chromosome in CHROMOSOMES:
                source = sources[class_name][chromosome]
                fixed_split = (
                    "validation" if chromosome == FIXED_SPLIT_CHROMOSOME["validation"]
                    else "test" if chromosome == FIXED_SPLIT_CHROMOSOME["test"]
                    else "train"
                )
                bed_handle = handles[
                    f"{class_name}_bed" if fixed_split == "train"
                    else f"{fixed_split}_{class_name}_bed"
                ]
                fasta_handle = handles[
                    f"{class_name}_fa" if fixed_split == "train"
                    else f"{fixed_split}_{class_name}_fa"
                ]
                for index in selections[class_name][chromosome]:
                    bed_record = source.bed_records[index]
                    fasta_record = source.fasta_records[index]
                    bed_handle.write(bed_record)
                    fasta_handle.write(fasta_record.text)
                    for destination in record_destinations(chromosome):
                        handles[f"{destination}_{class_name}_bed"].write(bed_record)
                        handles[f"{destination}_{class_name}_fa"].write(fasta_record.text)
                    bed_chrom, bed_start, bed_end = bed_coordinates(
                        bed_record, source.bed_path, index
                    )
                    manifest.writerow(
                        (class_name, fixed_split, f"chr{chromosome}",
                         CHROMOSOME_FOLD[chromosome], index + 1, bed_chrom,
                         bed_start, bed_end, fasta_record.identifier)
                    )

        training_count = sum(quotas[chrom] for chrom in TRAINING_CHROMOSOMES)
        validation_count = quotas[FIXED_SPLIT_CHROMOSOME["validation"]]
        test_count = quotas[FIXED_SPLIT_CHROMOSOME["test"]]
        summary = {
            "selected_benchmark_per_class": target,
            "selected_benchmark_total_records": target * 2,
            "training_per_class": training_count,
            "training_total_records": training_count * 2,
            "validation_per_class": validation_count,
            "validation_total_records": validation_count * 2,
            "test_per_class": test_count,
            "test_total_records": test_count * 2,
            "seed": seed,
            "excluded_chromosomes_from_training": ["chr1", "chr2", "chrY"],
            "fixed_split_chromosomes": {"validation": "chr1", "test": "chr2"},
            "quota_method": "Hamilton largest remainder across chr1-chr22 and chrX before fixed-split routing",
            "bed_fasta_pairing": "record position within each chromosome source pair",
            "folds": {
                f"fold{fold}": sum(quotas[chrom] for chrom in chromosomes)
                for fold, chromosomes in FOLD_CHROMOSOMES.items()
            },
            "outers": {
                f"outer{outer}": {
                    "train_per_class": sum(
                        quotas[chrom]
                        for fold in included_folds
                        for chrom in FOLD_CHROMOSOMES[fold]
                        if chrom != VALIDATION_CHROMOSOME[outer]
                    ),
                    "validation_per_class": quotas[VALIDATION_CHROMOSOME[outer]],
                    "held_out_fold": outer,
                }
                for outer, included_folds in OUTER_FOLDS.items()
            },
            "chromosomes": {
                f"chr{chromosome}": {
                    "positive_available": len(sources["pos"][chromosome].bed_records),
                    "negative_available": len(sources["neg"][chromosome].bed_records),
                    "positive_selected": quotas[chromosome],
                    "negative_selected": quotas[chromosome],
                }
                for chromosome in CHROMOSOMES
            },
        }
        handles["summary"].write(json.dumps(summary, indent=2, sort_keys=True) + "\n")

        for handle in handles.values():
            handle.flush()
            os.fsync(handle.fileno())
            handle.close()
        for name, path in paths.items():
            os.replace(temporary_paths[name], path)
    except BaseException:
        for handle in handles.values():
            if not handle.closed:
                handle.close()
        for temporary in temporary_paths.values():
            temporary.unlink(missing_ok=True)
        raise


def main() -> None:
    args = parse_args()
    chroms_dir = args.chroms_dir.resolve()
    output_dir = args.output_dir.resolve()
    folds_dir = args.folds_dir.resolve()
    sources = {
        class_name: {
            chromosome: load_source(chroms_dir, chromosome, class_name)
            for chromosome in CHROMOSOMES
        }
        for class_name in ("pos", "neg")
    }

    positive_counts = {
        chromosome: len(sources["pos"][chromosome].bed_records)
        for chromosome in CHROMOSOMES
    }
    quotas = hamilton_quotas(positive_counts, args.target_count)
    selections: dict[str, dict[str, tuple[int, ...]]] = {"pos": {}, "neg": {}}
    for class_offset, class_name in enumerate(("pos", "neg")):
        for chromosome_index, chromosome in enumerate(CHROMOSOMES):
            pool_size = len(sources[class_name][chromosome].bed_records)
            quota = quotas[chromosome]
            derived_seed = args.seed + class_offset * 100_000 + chromosome_index
            selections[class_name][chromosome] = sample_indices(
                pool_size, quota, derived_seed
            )

    paths = output_paths(output_dir, folds_dir)
    write_outputs(
        paths, sources, selections, quotas, args.target_count, args.seed, args.overwrite
    )
    training_count = sum(quotas[chromosome] for chromosome in TRAINING_CHROMOSOMES)
    print(
        f"Selected benchmark: {args.target_count:,} positives and "
        f"{args.target_count:,} negatives"
    )
    print(f"Fixed training pool: {training_count:,} per class")
    print(f"Validation (chr1): {quotas['1']:,} per class")
    print(f"Test (chr2): {quotas['2']:,} per class")
    for chromosome in CHROMOSOMES:
        print(f"chr{chromosome}: {quotas[chromosome]:,} positive; "
              f"{quotas[chromosome]:,} negative")
    for name, path in paths.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
