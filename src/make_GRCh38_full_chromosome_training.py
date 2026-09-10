#!/usr/bin/env python3
"""Write persistent 100% EVN and SVN chromosome-training pools.

Both pools exclude chr1 (validation) and chr2 (test), and concatenate training
chromosomes in the fixed order chr3--chr22, chrX. EVN records come from the
already balanced/downsampled chromosome shards. For SVN, each chromosome uses
the smaller positive/negative count; excess positives are selected with a fixed
seed while preserving BED/FASTA pairing and original record order. Both
synthetic and real genomic negative pools are supported and written to separate
``training/<negative-kind>/`` branches.
"""

from __future__ import annotations

import argparse
import random
from dataclasses import dataclass
from pathlib import Path


CHROMS = tuple(str(i) for i in range(3, 23)) + ("X",)


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    text: str


def parse_fasta(path: Path) -> tuple[FastaRecord, ...]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty FASTA: {path}")
    records: list[FastaRecord] = []
    header: str | None = None
    lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header, "".join(lines)))
                header = line[1:].strip()
                lines = [line]
            elif header is not None:
                lines.append(line)
    if header is not None:
        text = "".join(lines)
        records.append(FastaRecord(header, text if text.endswith("\n") else text + "\n"))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    if len({record.identifier for record in records}) != len(records):
        raise ValueError(f"Duplicate FASTA identifiers in {path}")
    return tuple(records)


def parse_bed(path: Path) -> tuple[str, ...]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty BED: {path}")
    with path.open(encoding="utf-8") as handle:
        return tuple(line if line.endswith("\n") else line + "\n" for line in handle)


def paired(path_bed: Path, path_fa: Path) -> tuple[tuple[str, ...], tuple[FastaRecord, ...]]:
    beds, fastas = parse_bed(path_bed), parse_fasta(path_fa)
    if len(beds) != len(fastas):
        raise ValueError(f"BED/FASTA mismatch: {path_bed} ({len(beds)}) vs {path_fa} ({len(fastas)})")
    return beds, fastas


def write_records(path: Path, records: list[str], overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(f"Refusing to overwrite {path}; pass --overwrite")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(records)
    temporary.replace(path)


def negative_paths(directory: Path, stem: str, negative_kind: str) -> tuple[Path, Path]:
    prefix = "synthetic_neg1x" if negative_kind == "synthetic" else "real_neg1x"
    return directory / f"{prefix}_{stem}.bed", directory / f"{prefix}_{stem}.fa"


def build_evn(repo: Path, negative_kind: str) -> dict[str, list[str]]:
    root = repo / f"data/processed/evn/GRCh38/training/{negative_kind}/chromosomes"
    output = {"pos_bed": [], "pos_fa": [], "neg_bed": [], "neg_fa": []}
    for chrom in CHROMS:
        directory = root / f"chrom{chrom}"
        for cls in ("pos", "neg"):
            beds, fastas = paired(directory / f"train_{cls}.bed", directory / f"train_{cls}.fa")
            output[f"{cls}_bed"].extend(beds)
            output[f"{cls}_fa"].extend(record.text for record in fastas)
    return output


def build_svn(repo: Path, seed: int, overwrite: bool, negative_kind: str) -> dict[str, list[str]]:
    root = repo / "data/processed/svn/GRCh38/chroms"
    output = {"pos_bed": [], "pos_fa": [], "neg_bed": [], "neg_fa": []}
    for offset, chrom in enumerate(CHROMS):
        directory = root / f"chrom{chrom}"
        stem = f"chrom{chrom}_GRCh38_Silencer_ELS"
        pos_bed, pos_fa = paired(directory / f"{stem}.bed", directory / f"{stem}.fa")
        neg_bed_path, neg_fa_path = negative_paths(directory, stem, negative_kind)
        neg_bed, neg_fa = paired(neg_bed_path, neg_fa_path)
        quota = min(len(pos_bed), len(neg_bed))
        pos_indices = sorted(random.Random(seed + offset).sample(range(len(pos_bed)), quota))
        neg_indices = sorted(random.Random(seed + 100_000 + offset).sample(range(len(neg_bed)), quota))
        selected = {
            "pos_bed": [pos_bed[i] for i in pos_indices],
            "pos_fa": [pos_fa[i].text for i in pos_indices],
            "neg_bed": [neg_bed[i] for i in neg_indices],
            "neg_fa": [neg_fa[i].text for i in neg_indices],
        }
        shard_dir = repo / f"data/processed/svn/GRCh38/training/{negative_kind}/chromosomes/chrom{chrom}"
        for cls in ("pos", "neg"):
            write_records(shard_dir / f"train_{cls}.bed", selected[f"{cls}_bed"], overwrite)
            write_records(shard_dir / f"train_{cls}.fa", selected[f"{cls}_fa"], overwrite)
        for key in output:
            output[key].extend(selected[key])

    # Fixed evaluation sets use the same balancing rule and never enter training.
    for split_offset, (split, chrom) in enumerate((("validation", "1"), ("test", "2"))):
        directory = root / f"chrom{chrom}"
        stem = f"chrom{chrom}_GRCh38_Silencer_ELS"
        pos_bed, pos_fa = paired(directory / f"{stem}.bed", directory / f"{stem}.fa")
        neg_bed_path, neg_fa_path = negative_paths(directory, stem, negative_kind)
        neg_bed, neg_fa = paired(neg_bed_path, neg_fa_path)
        quota = min(len(pos_bed), len(neg_bed))
        pos_indices = sorted(random.Random(seed + 200_000 + split_offset).sample(range(len(pos_bed)), quota))
        neg_indices = sorted(random.Random(seed + 300_000 + split_offset).sample(range(len(neg_bed)), quota))
        destination = repo / f"data/processed/svn/GRCh38/training/{negative_kind}"
        write_records(destination / f"{split}_pos.bed", [pos_bed[i] for i in pos_indices], overwrite)
        write_records(destination / f"{split}_pos.fa", [pos_fa[i].text for i in pos_indices], overwrite)
        write_records(destination / f"{split}_neg.bed", [neg_bed[i] for i in neg_indices], overwrite)
        write_records(destination / f"{split}_neg.fa", [neg_fa[i].text for i in neg_indices], overwrite)
    return output


def main() -> None:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--negative-kind",
        choices=("synthetic", "real"),
        default="synthetic",
        help="Negative source and output namespace (default: %(default)s).",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    for dataset, records in (
        ("evn", build_evn(repo, args.negative_kind)),
        ("svn", build_svn(repo, args.seed, args.overwrite, args.negative_kind)),
    ):
        destination = repo / f"data/processed/{dataset}/GRCh38/training/{args.negative_kind}"
        for cls in ("pos", "neg"):
            write_records(destination / f"train_100_{cls}.bed", records[f"{cls}_bed"], args.overwrite)
            write_records(destination / f"train_100_{cls}.fa", records[f"{cls}_fa"], args.overwrite)
        n_pos = len(records["pos_bed"])
        n_neg = len(records["neg_bed"])
        if n_pos == 0 or n_pos != n_neg:
            raise AssertionError(f"{dataset}: unbalanced full set ({n_pos} pos, {n_neg} neg)")
        print(f"{dataset}: wrote {n_pos:,} positives and {n_neg:,} negatives")


if __name__ == "__main__":
    main()
