#!/usr/bin/env python3
"""Build chromosome-resolved EVN/SVN positive inputs for negative generation.

SVN is retained in full (chr1--chr22 and chrX). EVN is downsampled to the
same total number of loci using Hamilton allocation so that its original
chromosome proportions are preserved as closely as integer counts allow.
No negative files are read by this stage.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import shutil
import tempfile
from pathlib import Path
from typing import Iterator


CHROMOSOMES = tuple(str(i) for i in range(1, 23)) + ("X",)


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/processed/consolidation/GRCh38"),
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def input_path(repo: Path, task: str, chromosome: str) -> Path:
    directory = repo / f"data/processed/{task}/GRCh38/chroms/chrom{chromosome}"
    suffix = "ELS" if task == "evn" else "Silencer_ELS"
    return directory / f"chrom{chromosome}_GRCh38_{suffix}.bed"


def bed_records(path: Path, expected_chromosome: str) -> Iterator[str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty BED input: {path}")
    expected = f"chr{expected_chromosome}"
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(f"{path}:{line_number}: blank BED record")
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
                raise ValueError(f"{path}:{line_number}: invalid BED coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")
            yield line if line.endswith("\n") else line + "\n"


def count_records(path: Path, chromosome: str) -> int:
    return sum(1 for _ in bed_records(path, chromosome))


def hamilton_quotas(counts: dict[str, int], target: int) -> dict[str, int]:
    """Allocate *target* records while minimizing chromosome proportion error."""
    total = sum(counts.values())
    if target < 0 or target > total:
        raise ValueError(f"Cannot allocate {target:,} records from {total:,}")
    exact = {chrom: target * counts[chrom] / total for chrom in CHROMOSOMES}
    quotas = {chrom: int(exact[chrom]) for chrom in CHROMOSOMES}
    remainder = target - sum(quotas.values())
    order = {chrom: index for index, chrom in enumerate(CHROMOSOMES)}
    ranked = sorted(
        CHROMOSOMES,
        key=lambda chrom: (-(exact[chrom] - quotas[chrom]), order[chrom]),
    )
    for chromosome in ranked[:remainder]:
        quotas[chromosome] += 1
    if sum(quotas.values()) != target:
        raise AssertionError("Hamilton quotas do not sum to the requested target")
    return quotas


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


def select_records(
    source: Path,
    destination: Path,
    chromosome: str,
    selected_indices: set[int],
) -> None:
    temporary = temporary_path(destination)
    written = 0
    try:
        with temporary.open("w", encoding="utf-8", newline="") as output:
            for index, record in enumerate(bed_records(source, chromosome)):
                if index in selected_indices:
                    output.write(record)
                    written += 1
            output.flush()
            os.fsync(output.fileno())
        if written != len(selected_indices):
            raise AssertionError(
                f"{source}: selected {len(selected_indices):,} records but wrote {written:,}"
            )
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def copy_bed(source: Path, destination: Path) -> None:
    temporary = temporary_path(destination)
    try:
        shutil.copyfile(source, temporary)
        os.replace(temporary, destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def concatenate_beds(first: Path, second: Path, destination: Path) -> None:
    temporary = temporary_path(destination)
    try:
        with temporary.open("wb") as output:
            for source in (first, second):
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
    output_root = args.output_dir
    if not output_root.is_absolute():
        output_root = repo / output_root

    evn_sources = {chrom: input_path(repo, "evn", chrom) for chrom in CHROMOSOMES}
    svn_sources = {chrom: input_path(repo, "svn", chrom) for chrom in CHROMOSOMES}
    evn_counts = {
        chrom: count_records(evn_sources[chrom], chrom) for chrom in CHROMOSOMES
    }
    svn_counts = {
        chrom: count_records(svn_sources[chrom], chrom) for chrom in CHROMOSOMES
    }
    target = sum(svn_counts.values())
    quotas = hamilton_quotas(evn_counts, target)

    destinations: list[Path] = []
    for chromosome in CHROMOSOMES:
        directory = output_root / "chromosomes" / f"chrom{chromosome}"
        destinations.extend(
            directory / filename
            for filename in ("evn_pos.bed", "svn_pos.bed", "consolidation.bed")
        )
    destinations.extend((output_root / "manifest.tsv", output_root / "summary.json"))
    existing = [path for path in destinations if path.exists()]
    if existing and not args.overwrite:
        preview = ", ".join(str(path) for path in existing[:3])
        raise FileExistsError(f"Outputs already exist ({preview}); pass --overwrite")

    manifest_rows: list[tuple[object, ...]] = []
    for index, chromosome in enumerate(CHROMOSOMES):
        directory = output_root / "chromosomes" / f"chrom{chromosome}"
        evn_output = directory / "evn_pos.bed"
        svn_output = directory / "svn_pos.bed"
        combined_output = directory / "consolidation.bed"
        chosen = set(
            random.Random(args.seed + index).sample(
                range(evn_counts[chromosome]), quotas[chromosome]
            )
        )
        select_records(evn_sources[chromosome], evn_output, chromosome, chosen)
        copy_bed(svn_sources[chromosome], svn_output)
        concatenate_beds(evn_output, svn_output, combined_output)

        observed_evn = count_records(evn_output, chromosome)
        observed_svn = count_records(svn_output, chromosome)
        observed_combined = count_records(combined_output, chromosome)
        if observed_evn != quotas[chromosome]:
            raise AssertionError(f"chr{chromosome}: EVN output count mismatch")
        if observed_svn != svn_counts[chromosome]:
            raise AssertionError(f"chr{chromosome}: SVN output count mismatch")
        if observed_combined != observed_evn + observed_svn:
            raise AssertionError(f"chr{chromosome}: consolidation count mismatch")
        manifest_rows.append(
            (
                f"chr{chromosome}",
                evn_counts[chromosome],
                observed_evn,
                observed_svn,
                observed_combined,
                evn_sources[chromosome].relative_to(repo),
                svn_sources[chromosome].relative_to(repo),
            )
        )

    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "manifest.tsv"
    manifest_tmp = temporary_path(manifest_path)
    summary_path = output_root / "summary.json"
    summary_tmp = temporary_path(summary_path)
    try:
        with manifest_tmp.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                (
                    "chromosome",
                    "evn_available",
                    "evn_selected",
                    "svn_selected",
                    "consolidation_records",
                    "evn_source",
                    "svn_source",
                )
            )
            writer.writerows(manifest_rows)
            handle.flush()
            os.fsync(handle.fileno())
        summary = {
            "seed": args.seed,
            "chromosomes": [f"chr{chrom}" for chrom in CHROMOSOMES],
            "excluded_chromosomes": ["chrY"],
            "evn_available": sum(evn_counts.values()),
            "evn_selected": sum(quotas.values()),
            "svn_selected": target,
            "consolidation_records": 2 * target,
            "evn_sampling": "Hamilton chromosome-proportional allocation followed by seeded sampling",
        }
        with summary_tmp.open("w", encoding="utf-8", newline="") as handle:
            json.dump(summary, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(manifest_tmp, manifest_path)
        os.replace(summary_tmp, summary_path)
    except BaseException:
        manifest_tmp.unlink(missing_ok=True)
        summary_tmp.unlink(missing_ok=True)
        raise

    print(f"EVN available (chr1-22/X): {sum(evn_counts.values()):,}")
    print(f"EVN selected:              {sum(quotas.values()):,}")
    print(f"SVN selected:              {target:,}")
    print(f"Consolidated records:      {2 * target:,}")


if __name__ == "__main__":
    main()
