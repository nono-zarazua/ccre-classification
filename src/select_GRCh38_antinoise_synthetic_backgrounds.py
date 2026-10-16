#!/usr/bin/env python3
"""Select one globally unique mix0 background per EVN/SVN foreground.

mix0 retains the foreground header and appends ``_N<index>_H<rank>``.  N is
the zero-based index in a chromosome's consolidation FASTA; H is the candidate
rank.  EVN records precede SVN records in every consolidation.

Synthetic sequences have no genomic locus, so leakage controls operate on
sequence identity.  A sequence and its reverse complement share one canonical
identity.  Candidates identical to any foreground are rejected, then a maximum
bipartite matching assigns unique candidate identities to as many foregrounds
as possible.  Positives are emitted only with an assigned negative.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import os
import random
import re
import shutil
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO


CHROMOSOMES = tuple(str(i) for i in range(1, 23)) + ("X",)
TRAIN_CHROMOSOMES = tuple(c for c in CHROMOSOMES if c not in ("1", "2"))
SYNTHETIC_HEADER = re.compile(r"^(?P<foreground>.+?)\s+_N(?P<index>\d+)_H(?P<rank>\d+)$")
COMPLEMENT = str.maketrans("ACGT", "TGCA")


@dataclass(frozen=True)
class BedRecord:
    chrom: str
    start: int
    end: int
    identifier: str
    text: str


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def text(self) -> str:
        return f">{self.header}\n{self.sequence}\n"


@dataclass(frozen=True)
class CandidateRef:
    rank: int
    identity: bytes


@dataclass(frozen=True)
class ChromosomeInfo:
    foreground_count: int
    evn_count: int
    svn_count: int


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--chromosomes", nargs="+", default=list(CHROMOSOMES),
        help="Chromosomes without chr prefix (default: chr1--chr22 and chrX)",
    )
    parser.add_argument(
        "--output-base", type=Path, default=Path("data/processed"),
        help="Output base, primarily useful for isolated tests",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def bed_records(path: Path) -> Iterator[BedRecord]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty BED: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(f"{path}:{line_number}: blank BED record")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least three columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")
            yield BedRecord(
                fields[0], start, end, fields[3] if len(fields) > 3 else "",
                line if line.endswith("\n") else line + "\n",
            )


def fasta_records(path: Path) -> Iterator[FastaRecord]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty FASTA: {path}")
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if header is not None:
                    if not sequence:
                        raise ValueError(f"{path}: empty sequence for {header}")
                    yield FastaRecord(header, "".join(sequence).upper())
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                sequence = []
            elif header is None:
                if line.strip():
                    raise ValueError(f"{path}:{line_number}: sequence before header")
            else:
                sequence.append(line.strip())
    if header is not None:
        if not sequence:
            raise ValueError(f"{path}: empty sequence for {header}")
        yield FastaRecord(header, "".join(sequence).upper())


def canonical_identity(sequence: str) -> bytes:
    reverse_complement = sequence.translate(COMPLEMENT)[::-1]
    canonical = sequence if sequence <= reverse_complement else reverse_complement
    return hashlib.sha256(canonical.encode("ascii")).digest()


def nucleotide_counts(sequence: str) -> tuple[int, int, int, int]:
    counts = Counter(sequence)
    return tuple(counts[base] for base in "ACGT")


def split_for(chromosome: str) -> str:
    return "validation" if chromosome == "1" else "test" if chromosome == "2" else "train"


def chromosome_dir(repo: Path, chromosome: str) -> Path:
    return repo / f"data/processed/consolidation/GRCh38/chromosomes/chrom{chromosome}"


def load_foregrounds(directory: Path) -> tuple[list[BedRecord], list[FastaRecord], ChromosomeInfo]:
    evn = list(bed_records(directory / "evn_pos.bed"))
    svn = list(bed_records(directory / "svn_pos.bed"))
    consolidated = list(bed_records(directory / "consolidation.bed"))
    fasta = list(fasta_records(directory / "consolidation.fa"))
    if consolidated != evn + svn:
        raise ValueError(f"{directory}: consolidation.bed is not EVN followed by SVN")
    if len(fasta) != len(consolidated):
        raise ValueError(f"{directory}: consolidation BED/FASTA count mismatch")
    for index, (bed, record) in enumerate(zip(consolidated, fasta)):
        header_fields = record.header.split()
        if len(header_fields) < 3:
            raise ValueError(f"{directory}: malformed foreground header at index {index}")
        header_locus = (header_fields[0], int(header_fields[1]), int(header_fields[2]))
        # longext_many writes a one-based start while retaining the BED end.
        # The BED remains the genomic source of truth.
        expected_header_locus = (bed.chrom, bed.start + 1, bed.end)
        if header_locus != expected_header_locus:
            raise ValueError(f"{directory}: foreground BED/FASTA mismatch at index {index}")
        if len(record.sequence) != bed.end - bed.start:
            raise ValueError(f"{directory}: foreground length mismatch at index {index}")
    # Explicitly validate the task boundary requested by the benchmark design.
    if consolidated[len(evn) - 1] != evn[-1] or consolidated[len(evn)] != svn[0]:
        raise AssertionError(f"{directory}: EVN/SVN boundary validation failed")
    return consolidated, fasta, ChromosomeInfo(len(consolidated), len(evn), len(svn))


def candidate_validation(
    record: FastaRecord,
    foregrounds: list[FastaRecord],
    positive_identities: set[bytes],
) -> tuple[int | None, int | None, bytes | None, str]:
    match = SYNTHETIC_HEADER.match(record.header)
    if match is None:
        return None, None, None, "malformed_header"
    index, rank = int(match.group("index")), int(match.group("rank"))
    if index < 0 or index >= len(foregrounds):
        return index, rank, None, "foreground_index_out_of_range"
    foreground = foregrounds[index]
    if match.group("foreground") != foreground.header:
        return index, rank, None, "foreground_header_mismatch"
    if len(record.sequence) != len(foreground.sequence):
        return index, rank, None, "length_mismatch"
    if set(record.sequence) - set("ACGT"):
        return index, rank, None, "non_acgt_sequence"
    if nucleotide_counts(record.sequence) != nucleotide_counts(foreground.sequence):
        return index, rank, None, "mononucleotide_composition_mismatch"
    identity = canonical_identity(record.sequence)
    if identity in positive_identities:
        return index, rank, identity, "positive_sequence_collision"
    return index, rank, identity, "valid"


def maximum_matching(
    options: dict[tuple[str, int], list[CandidateRef]], seed: int
) -> dict[tuple[str, int], CandidateRef]:
    keys = list(options)
    random.Random(seed).shuffle(keys)
    keys.sort(key=lambda key: len(options[key]))
    owner: dict[bytes, tuple[str, int]] = {}
    chosen: dict[tuple[str, int], CandidateRef] = {}
    sys.setrecursionlimit(max(10_000, len(keys) + 100))

    def augment(key: tuple[str, int], visited: set[bytes]) -> bool:
        for candidate in sorted(options[key], key=lambda item: item.rank):
            if candidate.identity in visited:
                continue
            visited.add(candidate.identity)
            previous = owner.get(candidate.identity)
            if previous is None or augment(previous, visited):
                owner[candidate.identity] = key
                chosen[key] = candidate
                return True
        return False

    for key in keys:
        augment(key, set())
    return chosen


def open_output(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", encoding="utf-8", newline="")


def write_pair(
    handles: dict[str, TextIO], prefix: str, bed: BedRecord,
    positive: FastaRecord, negative: FastaRecord,
) -> None:
    handles[f"{prefix}_pos_bed"].write(bed.text)
    handles[f"{prefix}_pos_fa"].write(positive.text)
    handles[f"{prefix}_neg_fa"].write(negative.text)


def main() -> None:
    args = parse_args()
    repo = args.repo.resolve()
    chromosomes = tuple(args.chromosomes)
    if len(chromosomes) != len(set(chromosomes)) or set(chromosomes) - set(CHROMOSOMES):
        raise ValueError(f"Invalid or repeated chromosomes: {chromosomes}")
    output_base = args.output_base if args.output_base.is_absolute() else repo / args.output_base

    infos: dict[str, ChromosomeInfo] = {}
    positive_identities: set[bytes] = set()
    for chromosome in chromosomes:
        _, foregrounds, info = load_foregrounds(chromosome_dir(repo, chromosome))
        infos[chromosome] = info
        positive_identities.update(canonical_identity(record.sequence) for record in foregrounds)

    options: dict[tuple[str, int], list[CandidateRef]] = {}
    candidate_status_counts: Counter[str] = Counter()
    for chromosome in chromosomes:
        directory = chromosome_dir(repo, chromosome)
        _, foregrounds, _ = load_foregrounds(directory)
        grouped: dict[int, list[CandidateRef]] = defaultdict(list)
        ranks: dict[int, set[int]] = defaultdict(set)
        for record in fasta_records(directory / "synthetic_antinoise_candidates_2x.fa"):
            index, rank, identity, status = candidate_validation(record, foregrounds, positive_identities)
            candidate_status_counts[status] += 1
            if status != "valid":
                continue
            assert index is not None and rank is not None and identity is not None
            if rank in ranks[index]:
                candidate_status_counts["duplicate_candidate_rank"] += 1
                continue
            ranks[index].add(rank)
            grouped[index].append(CandidateRef(rank, identity))
        for index in range(infos[chromosome].foreground_count):
            options[(chromosome, index)] = grouped.get(index, [])

    chosen = maximum_matching(options, args.seed)
    task_roots = {task: output_base / task / "GRCh38/synthetic" for task in ("evn", "svn")}
    audit_root = output_base / "consolidation/GRCh38/antinoise_synthetic_selection"
    destinations = [*task_roots.values(), audit_root]
    existing = [path for path in destinations if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Output exists: {existing[0]}; pass --overwrite")
    output_base.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="antinoise-synthetic-selection.", dir=output_base))
    staged_roots = {task: staging / task / "GRCh38/synthetic" for task in ("evn", "svn")}
    staged_audit = staging / "consolidation/GRCh38/antinoise_synthetic_selection"
    handles: dict[str, TextIO] = {}
    counts: dict[str, Counter[str]] = {task: Counter() for task in ("evn", "svn")}
    try:
        for task, root in staged_roots.items():
            for split in ("validation", "test", "train_100"):
                handles[f"{task}_{split}_pos_bed"] = open_output(root / f"{split}_pos.bed")
                handles[f"{task}_{split}_pos_fa"] = open_output(root / f"{split}_pos.fa")
                handles[f"{task}_{split}_neg_fa"] = open_output(root / f"{split}_neg.fa")
        provenance_handle = open_output(staged_audit / "provenance.tsv")
        provenance = csv.writer(provenance_handle, delimiter="\t", lineterminator="\n")
        provenance.writerow((
            "chromosome", "split", "task", "foreground_index_zero_based",
            "foreground_row_one_based", "foreground_id", "candidate_count",
            "selected_rank", "status",
        ))

        for chromosome in chromosomes:
            directory = chromosome_dir(repo, chromosome)
            beds, foregrounds, info = load_foregrounds(directory)
            selected_records: dict[int, FastaRecord] = {}
            selected_refs = {index: chosen[(chromosome, index)] for index in range(info.foreground_count) if (chromosome, index) in chosen}
            for record in fasta_records(directory / "synthetic_antinoise_candidates_2x.fa"):
                index, rank, identity, status = candidate_validation(
                    record, foregrounds, positive_identities
                )
                if status != "valid" or index is None or rank is None:
                    continue
                selected = selected_refs.get(index)
                if selected is not None and selected.rank == rank and selected.identity == identity:
                    selected_records[index] = record
            split = split_for(chromosome)
            for index, (bed, positive) in enumerate(zip(beds, foregrounds)):
                task = "evn" if index < info.evn_count else "svn"
                negative = selected_records.get(index)
                if negative is None:
                    status = "no_candidate_generated" if not options[(chromosome, index)] else "no_unique_candidate_assignment"
                    provenance.writerow((f"chr{chromosome}", split, task, index, index + 1, bed.identifier, len(options[(chromosome, index)]), "", status))
                    continue
                destination = "train_100" if split == "train" else split
                write_pair(handles, f"{task}_{destination}", bed, positive, negative)
                if split == "train":
                    chrom_root = staged_roots[task] / "chromosomes" / f"chrom{chromosome}"
                    prefix = f"{task}_chrom{chromosome}"
                    if f"{prefix}_pos_bed" not in handles:
                        handles[f"{prefix}_pos_bed"] = open_output(chrom_root / "train_pos.bed")
                        handles[f"{prefix}_pos_fa"] = open_output(chrom_root / "train_pos.fa")
                        handles[f"{prefix}_neg_fa"] = open_output(chrom_root / "train_neg.fa")
                    write_pair(handles, prefix, bed, positive, negative)
                counts[task][split] += 1
                provenance.writerow((f"chr{chromosome}", split, task, index, index + 1, bed.identifier, len(options[(chromosome, index)]), selected_refs[index].rank, "selected"))

        for handle in handles.values():
            handle.flush(); os.fsync(handle.fileno()); handle.close()
        provenance_handle.flush(); os.fsync(provenance_handle.fileno()); provenance_handle.close()
        with (staged_audit / "summary.json").open("w", encoding="utf-8") as handle:
            json.dump({
                "seed": args.seed,
                "chromosomes": [f"chr{chromosome}" for chromosome in chromosomes],
                "identity_rule": "SHA-256 of lexicographically canonical sequence/reverse-complement",
                "selection": "maximum bipartite matching; H1 preferred during deterministic augmentation",
                "candidate_status_counts": dict(sorted(candidate_status_counts.items())),
                "foregrounds": len(options),
                "selected_pairs": len(chosen),
                "unresolved_foregrounds": len(options) - len(chosen),
                "records_per_class": {task: dict(value) for task, value in counts.items()},
            }, handle, indent=2, sort_keys=True)
            handle.write("\n")

        for destination in destinations:
            source = staging / destination.relative_to(output_base)
            if destination.exists():
                shutil.rmtree(destination)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(source), str(destination))
    except BaseException:
        for handle in handles.values():
            if not handle.closed:
                handle.close()
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    for task in ("evn", "svn"):
        print(
            f"{task.upper()}: validation={counts[task]['validation']:,}; "
            f"test={counts[task]['test']:,}; training={counts[task]['train']:,}"
        )
    print(f"Selected {len(chosen):,} of {len(options):,} foreground/background pairs")
    print(f"Audit: {audit_root}")


if __name__ == "__main__":
    main()
