#!/usr/bin/env python3
"""Build reproducible, strictly nested subsets from the fixed master pool.

Positive and negative BED/FASTA pairs are shuffled independently exactly once.
Every learning-curve dataset is then a prefix of that one fixed order, which
guarantees 10% is contained in 20%, 20% in 40%, and so on through 100%.
"""

from __future__ import annotations
import argparse
import json
import os
import random
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

PERCENTAGES = (10, 20, 40, 60, 80, 100)
CLASSES = ("pos", "neg")
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]

@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    text: str

@dataclass(frozen=True)
class PairedRecord:
    bed: str
    fasta: FastaRecord

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create nested GRCh38 learning-curve BED/FASTA datasets.")
    parser.add_argument("--input-dir", type=Path, default=REPOSITORY_ROOT / "data/processed/evn/GRCh38/training", help="Directory containing train_100_{pos,neg}.{bed,fa}.")
    parser.add_argument("--output-dir", type=Path, default=REPOSITORY_ROOT / "data/processed/evn/GRCh38/learning_curves")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true", help="Replace existing outputs; otherwise existing files cause failure.")
    return parser.parse_args()

def read_bed(path: Path) -> list[str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"BED input is missing or empty: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        records = [line if line.endswith("\n") else line + "\n" for line in handle]
    if any(not record.strip() for record in records):
        raise ValueError(f"Blank BED record found in {path}")
    return records

def read_fasta(path: Path) -> list[FastaRecord]:
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
                        raise ValueError(f"{path}: {identifier!r} has no sequence")
                    text = "".join(lines)
                    records.append(FastaRecord(identifier, text if text.endswith("\n") else text + "\n"))
                identifier = line[1:].strip().split(maxsplit=1)[0]
                if not identifier:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                if identifier in seen:
                    raise ValueError(f"{path}:{line_number}: duplicate ID {identifier!r}")
                seen.add(identifier)
                lines = [line]
            elif identifier is None:
                if line.strip():
                    raise ValueError(f"{path}:{line_number}: sequence before first header")
            else:
                lines.append(line)
    if identifier is not None:
        if len(lines) == 1:
            raise ValueError(f"{path}: {identifier!r} has no sequence")
        text = "".join(lines)
        records.append(FastaRecord(identifier, text if text.endswith("\n") else text + "\n"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records

def load_pairs(input_dir: Path, class_name: str) -> list[PairedRecord]:
    beds = read_bed(input_dir / f"train_100_{class_name}.bed")
    fastas = read_fasta(input_dir / f"train_100_{class_name}.fa")
    if len(beds) != len(fastas):
        raise ValueError(f"Master {class_name} BED/FASTA mismatch: {len(beds):,} versus {len(fastas):,}")
    return [PairedRecord(bed, fasta) for bed, fasta in zip(beds, fastas)]

def subset_sizes(total: int) -> dict[int, int]:
    """Use floor(total * percentage / 100), capped at total."""
    sizes = {percentage: min(total, total * percentage // 100) for percentage in PERCENTAGES}
    sizes[100] = total
    if list(sizes.values()) != sorted(sizes.values()):
        raise AssertionError(f"Non-monotonic subset sizes: {sizes}")
    return sizes

def output_paths(output_dir: Path) -> list[Path]:
    paths = [output_dir / f"train_{percentage}_{class_name}.{extension}" for percentage in PERCENTAGES for class_name in CLASSES for extension in ("bed", "fa")]
    return paths + [output_dir / "learning_curve_summary.json"]

def ensure_outputs_are_safe(paths: Sequence[Path], overwrite: bool) -> None:
    existing = [path for path in paths if path.exists()]
    if existing and not overwrite:
        listing = "\n".join(f"  - {path}" for path in existing)
        raise FileExistsError(f"Outputs already exist; use --overwrite to replace them:\n{listing}")

def write_atomic(path: Path, records: Sequence[PairedRecord], kind: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary: str | None = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False) as handle:
            temporary = handle.name
            for record in records:
                handle.write(record.bed if kind == "bed" else record.fasta.text)
        os.replace(temporary, path)
    except BaseException:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)
        raise

def validate_prefixes(prefixes: dict[int, list[PairedRecord]], class_name: str) -> None:
    counts = [len(prefixes[p]) for p in PERCENTAGES]
    if counts != sorted(counts):
        raise AssertionError(f"{class_name} counts are not monotonic: {counts}")
    for smaller, larger in zip(PERCENTAGES, PERCENTAGES[1:]):
        small_records = prefixes[smaller]
        if small_records != prefixes[larger][: len(small_records)]:
            raise AssertionError(f"{class_name}: {smaller}% is not an exact prefix of {larger}%")
    ids = [record.fasta.identifier for record in prefixes[100]]
    if len(ids) != len(set(ids)):
        raise AssertionError(f"Duplicate {class_name} FASTA IDs found after shuffling")

def chromosome_counts(records: Sequence[PairedRecord]) -> dict[str, int]:
    counts = Counter(record.bed.split("\t", 1)[0] for record in records)
    return dict(sorted(counts.items()))


def validate_written(output_dir: Path, class_name: str, sizes: dict[int, int]) -> None:
    previous_beds: list[str] = []
    previous_ids: list[str] = []
    for percentage in PERCENTAGES:
        beds = read_bed(output_dir / f"train_{percentage}_{class_name}.bed")
        fastas = read_fasta(output_dir / f"train_{percentage}_{class_name}.fa")
        ids = [record.identifier for record in fastas]
        expected = sizes[percentage]
        if len(beds) != expected or len(ids) != expected:
            raise AssertionError(
                f"{class_name} {percentage}%: expected {expected:,} synchronized "
                f"records; observed {len(beds):,} BED and {len(ids):,} FASTA"
            )
        if previous_beds != beds[: len(previous_beds)]:
            raise AssertionError(f"{class_name} BED nesting failed at {percentage}%")
        if previous_ids != ids[: len(previous_ids)]:
            raise AssertionError(f"{class_name} FASTA nesting failed at {percentage}%")
        previous_beds, previous_ids = beds, ids


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    if input_dir == output_dir:
        raise ValueError("Input and output directories must differ")
    ensure_outputs_are_safe(output_paths(output_dir), args.overwrite)
    summary: dict[str, object] = {
        "seed": args.seed,
        "percentages": list(PERCENTAGES),
        "rounding": "floor(total_records * percentage / 100); 100% is exactly total_records",
        "strategy": "one independent shuffle per class; every subset is an exact prefix",
        "classes": {},
    }
    class_summaries: dict[str, object] = {}
    for class_offset, class_name in enumerate(CLASSES):
        records = load_pairs(input_dir, class_name)
        random.Random(args.seed + class_offset).shuffle(records)
        sizes = subset_sizes(len(records))
        prefixes = {percentage: records[: sizes[percentage]] for percentage in PERCENTAGES}
        validate_prefixes(prefixes, class_name)
        subset_summary: dict[str, object] = {}
        for percentage in PERCENTAGES:
            prefix = prefixes[percentage]
            for extension in ("bed", "fa"):
                write_atomic(output_dir / f"train_{percentage}_{class_name}.{extension}", prefix, extension)
            subset_summary[str(percentage)] = {"records": len(prefix), "chromosome_counts": chromosome_counts(prefix)}
        validate_written(output_dir, class_name, sizes)
        class_summaries[class_name] = subset_summary
        print(f"{class_name}: " + ", ".join(f"{p}%={sizes[p]:,}" for p in PERCENTAGES))
    summary["classes"] = class_summaries
    summary_path = output_dir / "learning_curve_summary.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=summary_path.parent, prefix=f".{summary_path.name}.", suffix=".tmp", delete=False) as handle:
        temporary = handle.name
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, summary_path)
    print("Nested BED/FASTA prefix validation: PASS")

if __name__ == "__main__":
    main()
