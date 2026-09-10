#!/usr/bin/env python3
"""Build a chromosome-proportional EVN benchmark for chromosome-wise training.

The total number of EVN positives is matched to the cleaned SVN positive set.
That target is allocated across chr1--chr22 and chrX according to the chromosome
distribution of the complete clean EVN collection. BED and FASTA records are
paired by position, sampled once per chromosome and class, and kept together.

chr1 is written as the fixed validation set, chr2 as the fixed test set, and
chr3--chr22 plus chrX as separate training shards. No percentage subsets,
folds, outer sets, or cumulative chromosome files are produced.
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
VALIDATION_CHROMOSOME = "1"
TEST_CHROMOSOME = "2"
TRAINING_CHROMOSOMES = tuple(
    chromosome
    for chromosome in CHROMOSOMES
    if chromosome not in (VALIDATION_CHROMOSOME, TEST_CHROMOSOME)
)


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    text: str


@dataclass(frozen=True)
class PairedSource:
    bed_path: Path
    fasta_path: Path
    bed_records: tuple[str, ...]
    fasta_records: tuple[FastaRecord, ...]


def parse_args() -> argparse.Namespace:
    repository = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "Match clean EVN to the cleaned SVN locus count while preserving "
            "EVN chromosome proportions, then write chromosome training shards."
        )
    )
    parser.add_argument(
        "--chroms-dir",
        type=Path,
        default=repository / "data/processed/evn/GRCh38/chroms",
        help="EVN chromosome directory containing synchronized positive/negative BED/FASTA files.",
    )
    parser.add_argument(
        "--svn-bed",
        type=Path,
        default=repository / "data/processed/svn/GRCh38/GRCh38_Silencer_ELS.trimmed.bed",
        help="Cleaned SVN BED whose record count defines the EVN target size.",
    )
    parser.add_argument(
        "--svn-chroms-dir",
        type=Path,
        default=repository / "data/processed/svn/GRCh38/chroms",
        help="SVN chromosome sources used to derive the balanced real-background target.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory. Defaults to data/processed/evn/GRCh38/training/"
            "<negative-kind>."
        ),
    )
    parser.add_argument(
        "--negative-kind",
        choices=("synthetic", "real"),
        default="synthetic",
        help="Select genNullSeqs synthetic negatives or antinoise real genomic backgrounds.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing outputs. Existing files otherwise cause a safe failure.",
    )
    return parser.parse_args()


def source_paths(
    chroms_dir: Path,
    chromosome: str,
    class_name: str,
    negative_kind: str,
) -> tuple[Path, Path]:
    directory = chroms_dir / f"chrom{chromosome}"
    positive_stem = f"chrom{chromosome}_GRCh38_ELS"
    if class_name == "pos":
        return directory / f"{positive_stem}.bed", directory / f"{positive_stem}.fa"

    prefix = "synthetic_neg1x" if negative_kind == "synthetic" else "real_neg1x"
    bed_path = directory / f"{prefix}_{positive_stem}.bed"
    fasta_path = directory / f"{prefix}_{positive_stem}.fa"
    # gkmSVM::genNullSeqs() drops the literal x from xfold=1 FASTA names in
    # the current environment, although it retains it in the requested BED name.
    if negative_kind == "synthetic" and not fasta_path.is_file():
        fasta_path = directory / f"synthetic_neg1_{positive_stem}.fa"
    return bed_path, fasta_path


def read_bed(path: Path, *, allow_empty: bool = False) -> tuple[str, ...]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing BED input: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        records = tuple(line if line.endswith("\n") else line + "\n" for line in handle)
    if not records and not allow_empty:
        raise ValueError(f"BED input is empty: {path}")
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
                        raise ValueError(f"{path}: {identifier!r} has no sequence")
                    text = "".join(lines)
                    records.append(FastaRecord(identifier, text if text.endswith("\n") else text + "\n"))
                identifier = line[1:].strip()
                if not identifier:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                if identifier in seen:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA ID {identifier!r}")
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
    return tuple(records)


def load_source(
    chroms_dir: Path,
    chromosome: str,
    class_name: str,
    negative_kind: str,
) -> PairedSource:
    bed_path, fasta_path = source_paths(
        chroms_dir, chromosome, class_name, negative_kind
    )
    beds = read_bed(bed_path)
    fastas = read_fasta(fasta_path)
    if len(beds) != len(fastas):
        raise ValueError(
            f"chr{chromosome} {class_name} BED/FASTA mismatch: "
            f"{len(beds):,} BED versus {len(fastas):,} FASTA records"
        )
    return PairedSource(bed_path, fasta_path, beds, fastas)


def hamilton_quotas(counts: dict[str, int], target: int) -> dict[str, int]:
    """Allocate an exact target using proportional largest-remainder rounding."""
    total = sum(counts.values())
    if target <= 0:
        raise ValueError("SVN target count must be positive")
    if target > total:
        raise ValueError(f"SVN target {target:,} exceeds {total:,} available EVN positives")
    exact = {chromosome: target * count / total for chromosome, count in counts.items()}
    quotas = {chromosome: int(value) for chromosome, value in exact.items()}
    remaining = target - sum(quotas.values())
    order = {chromosome: index for index, chromosome in enumerate(CHROMOSOMES)}
    ranked = sorted(
        counts,
        key=lambda chromosome: (
            -(exact[chromosome] - quotas[chromosome]), order[chromosome]
        ),
    )
    for chromosome in ranked[:remaining]:
        quotas[chromosome] += 1
    if sum(quotas.values()) != target:
        raise AssertionError("Hamilton allocation does not sum to the SVN target")
    return quotas


def sample_indices(size: int, quota: int, seed: int) -> tuple[int, ...]:
    if quota > size:
        raise ValueError(f"Requested {quota:,} records from a pool of {size:,}")
    # Sorted indices randomize membership while preserving source order in outputs.
    return tuple(sorted(random.Random(seed).sample(range(size), quota)))


def output_paths(output_dir: Path) -> dict[str, Path]:
    paths: dict[str, Path] = {
        "validation_pos_bed": output_dir / "validation_pos.bed",
        "validation_pos_fa": output_dir / "validation_pos.fa",
        "validation_neg_bed": output_dir / "validation_neg.bed",
        "validation_neg_fa": output_dir / "validation_neg.fa",
        "test_pos_bed": output_dir / "test_pos.bed",
        "test_pos_fa": output_dir / "test_pos.fa",
        "test_neg_bed": output_dir / "test_neg.bed",
        "test_neg_fa": output_dir / "test_neg.fa",
        "manifest": output_dir / "chromosome_training_manifest.tsv",
        "summary": output_dir / "chromosome_training_summary.json",
    }
    for chromosome in TRAINING_CHROMOSOMES:
        directory = output_dir / "chromosomes" / f"chrom{chromosome}"
        for class_name in ("pos", "neg"):
            paths[f"chrom{chromosome}_{class_name}_bed"] = directory / f"train_{class_name}.bed"
            paths[f"chrom{chromosome}_{class_name}_fa"] = directory / f"train_{class_name}.fa"
    return paths


def temporary_handle(path: Path) -> tuple[Path, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", dir=path.parent,
        prefix=f".{path.name}.", suffix=".tmp", delete=False,
    )
    return Path(handle.name), handle


def destination(chromosome: str, class_name: str, extension: str) -> str:
    if chromosome == VALIDATION_CHROMOSOME:
        return f"validation_{class_name}_{extension}"
    if chromosome == TEST_CHROMOSOME:
        return f"test_{class_name}_{extension}"
    return f"chrom{chromosome}_{class_name}_{extension}"


def write_outputs(
    paths: dict[str, Path],
    sources: dict[str, dict[str, PairedSource]],
    selections: dict[str, dict[str, tuple[int, ...]]],
    quotas: dict[str, int],
    target: int,
    seed: int,
    negative_kind: str,
    overwrite: bool,
) -> None:
    existing = [path for path in paths.values() if path.exists()]
    if existing and not overwrite:
        raise FileExistsError(
            "Refusing to overwrite existing output(s); pass --overwrite:\n  "
            + "\n  ".join(str(path) for path in existing)
        )

    temporary_paths: dict[str, Path] = {}
    handles: dict[str, object] = {}
    try:
        for name, path in paths.items():
            temporary_paths[name], handles[name] = temporary_handle(path)

        manifest = csv.writer(handles["manifest"], delimiter="\t", lineterminator="\n")
        manifest.writerow(
            ("class", "split", "source_chromosome", "source_record_1based", "fasta_id")
        )
        for class_name in ("pos", "neg"):
            for chromosome in CHROMOSOMES:
                source = sources[class_name][chromosome]
                split = (
                    "validation" if chromosome == VALIDATION_CHROMOSOME
                    else "test" if chromosome == TEST_CHROMOSOME
                    else f"train_chr{chromosome}"
                )
                bed_handle = handles[destination(chromosome, class_name, "bed")]
                fasta_handle = handles[destination(chromosome, class_name, "fa")]
                for index in selections[class_name][chromosome]:
                    bed_handle.write(source.bed_records[index])
                    fasta_handle.write(source.fasta_records[index].text)
                    manifest.writerow(
                        (class_name, split, f"chr{chromosome}", index + 1,
                         source.fasta_records[index].identifier)
                    )

        summary = {
            "strategy": "independent chromosome training shards; no percentage subsets",
            "target_source": (
                "sum of per-chromosome min(SVN positives, real negatives)"
                if negative_kind == "real"
                else "record count of the cleaned SVN trimmed BED"
            ),
            "selected_benchmark_per_class": target,
            "selected_benchmark_total_records": target * 2,
            "seed": seed,
            "negative_kind": negative_kind,
            "fixed_split_chromosomes": {"validation": "chr1", "test": "chr2"},
            "training_chromosome_order": [f"chr{chrom}" for chrom in TRAINING_CHROMOSOMES],
            "quota_method": "Hamilton largest remainder over clean EVN chr1-chr22 and chrX",
            "bed_fasta_pairing": "record position within each chromosome and class",
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


def validate_outputs(paths: dict[str, Path], quotas: dict[str, int]) -> None:
    for class_name in ("pos", "neg"):
        for chromosome in CHROMOSOMES:
            beds = read_bed(paths[destination(chromosome, class_name, "bed")])
            fastas = read_fasta(paths[destination(chromosome, class_name, "fa")])
            expected = quotas[chromosome]
            if len(beds) != expected or len(fastas) != expected:
                raise AssertionError(
                    f"chr{chromosome} {class_name}: expected {expected:,}; "
                    f"found {len(beds):,} BED and {len(fastas):,} FASTA records"
                )
            if len({record.identifier for record in fastas}) != len(fastas):
                raise AssertionError(f"Duplicate FASTA IDs in chr{chromosome} {class_name}")


def balanced_real_svn_target(svn_chroms_dir: Path) -> int:
    """Return the balanced SVN size supported by real negatives on every chromosome."""
    target = 0
    for chromosome in CHROMOSOMES:
        directory = svn_chroms_dir / f"chrom{chromosome}"
        stem = f"chrom{chromosome}_GRCh38_Silencer_ELS"
        positive_bed = read_bed(directory / f"{stem}.bed")
        positive_fasta = read_fasta(directory / f"{stem}.fa")
        negative_bed = read_bed(directory / f"real_neg1x_{stem}.bed")
        negative_fasta = read_fasta(directory / f"real_neg1x_{stem}.fa")
        if len(positive_bed) != len(positive_fasta):
            raise ValueError(f"chr{chromosome} SVN positive BED/FASTA mismatch")
        if len(negative_bed) != len(negative_fasta):
            raise ValueError(f"chr{chromosome} SVN real-negative BED/FASTA mismatch")
        target += min(len(positive_bed), len(negative_bed))
    return target


def main() -> None:
    args = parse_args()
    chroms_dir = args.chroms_dir.resolve()
    svn_bed = args.svn_bed.resolve()
    output_dir = (
        args.output_dir
        if args.output_dir is not None
        else Path(__file__).resolve().parents[1]
        / "data/processed/evn/GRCh38/training"
        / args.negative_kind
    ).resolve()
    if args.negative_kind == "real":
        target = balanced_real_svn_target(args.svn_chroms_dir.resolve())
    else:
        target = len(read_bed(svn_bed))

    sources = {
        class_name: {
            chromosome: load_source(
                chroms_dir, chromosome, class_name, args.negative_kind
            )
            for chromosome in CHROMOSOMES
        }
        for class_name in ("pos", "neg")
    }
    positive_counts = {
        chromosome: len(sources["pos"][chromosome].bed_records)
        for chromosome in CHROMOSOMES
    }
    quotas = hamilton_quotas(positive_counts, target)

    selections: dict[str, dict[str, tuple[int, ...]]] = {"pos": {}, "neg": {}}
    for class_offset, class_name in enumerate(("pos", "neg")):
        for chromosome_index, chromosome in enumerate(CHROMOSOMES):
            selections[class_name][chromosome] = sample_indices(
                len(sources[class_name][chromosome].bed_records),
                quotas[chromosome],
                args.seed + class_offset * 100_000 + chromosome_index,
            )

    paths = output_paths(output_dir)
    write_outputs(
        paths, sources, selections, quotas, target, args.seed,
        args.negative_kind, args.overwrite
    )
    validate_outputs(paths, quotas)

    print(f"Selected EVN benchmark: {target:,} positives and {target:,} negatives")
    print(f"Negative source: {args.negative_kind}")
    print(f"Validation chr1: {quotas['1']:,} per class")
    print(f"Test chr2: {quotas['2']:,} per class")
    for chromosome in TRAINING_CHROMOSOMES:
        print(f"Training chr{chromosome}: {quotas[chromosome]:,} per class")
    print("BED/FASTA synchronization validation: PASS")


if __name__ == "__main__":
    main()
