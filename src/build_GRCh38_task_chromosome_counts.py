#!/usr/bin/env python3
"""Build current EVN/SVN chromosome tables from consolidated selections."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from collections import Counter
from pathlib import Path


TASKS = ("evn", "svn")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("data/processed"),
        help="Processed-data root (default: data/processed)",
    )
    return parser.parse_args()


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty consolidation manifest: {path}")
    return rows


def selected_counts(path: Path) -> tuple[Counter[tuple[str, str]], Counter[tuple[str, str]]]:
    selected: Counter[tuple[str, str]] = Counter()
    total: Counter[tuple[str, str]] = Counter()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chromosome", "task", "status"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required columns in {path}")
        for row in reader:
            key = (row["task"], row["chromosome"])
            total[key] += 1
            if row["status"] == "selected":
                selected[key] += 1
    return selected, total


def split_for(chromosome: str) -> str:
    if chromosome == "chr1":
        return "validation"
    if chromosome == "chr2":
        return "test"
    return "train"


def atomic_write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = (
        "chromosome",
        "split",
        "clean_positive_loci",
        "benchmark_foregrounds",
        "real_pairs",
        "missing_real_pairs",
        "synthetic_pairs",
        "missing_synthetic_pairs",
    )
    handle = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    )
    temporary = Path(handle.name)
    try:
        with handle:
            writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> None:
    args = parse_args()
    consolidation = args.root / "consolidation" / "GRCh38"
    manifest = read_manifest(consolidation / "manifest.tsv")
    real_selected, real_total = selected_counts(
        consolidation / "antinoise_selection" / "provenance.tsv"
    )
    synthetic_selected, synthetic_total = selected_counts(
        consolidation / "antinoise_synthetic_selection" / "provenance.tsv"
    )

    for task in TASKS:
        rows: list[dict[str, object]] = []
        for source in manifest:
            chromosome = source["chromosome"]
            clean = int(source["evn_available"] if task == "evn" else source["svn_selected"])
            benchmark = int(source["evn_selected"] if task == "evn" else source["svn_selected"])
            key = (task, chromosome)
            if real_total[key] != benchmark:
                raise ValueError(
                    f"{task} {chromosome}: real provenance has {real_total[key]} "
                    f"rows; expected {benchmark}"
                )
            if synthetic_total[key] != benchmark:
                raise ValueError(
                    f"{task} {chromosome}: synthetic provenance has "
                    f"{synthetic_total[key]} rows; expected {benchmark}"
                )
            real = real_selected[key]
            synthetic = synthetic_selected[key]
            rows.append(
                {
                    "chromosome": chromosome,
                    "split": split_for(chromosome),
                    "clean_positive_loci": clean,
                    "benchmark_foregrounds": benchmark,
                    "real_pairs": real,
                    "missing_real_pairs": benchmark - real,
                    "synthetic_pairs": synthetic,
                    "missing_synthetic_pairs": benchmark - synthetic,
                }
            )

        output = args.root / task / "GRCh38" / "chroms" / "chromosome_counts.csv"
        atomic_write(output, rows)
        print(f"Wrote {len(rows)} chromosome rows to {output}")


if __name__ == "__main__":
    main()
