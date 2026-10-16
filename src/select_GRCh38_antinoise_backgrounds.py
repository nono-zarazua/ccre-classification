#!/usr/bin/env python3
"""Select one antinoise background per EVN/SVN foreground.

Antinoise labels a candidate ``peakN_R`` where N is the zero-based record
index in the chromosome's consolidated foreground FASTA and R is candidate
rank 1 or 2.  Candidate FASTA and BED records are synchronized by position,
but candidate output order is otherwise arbitrary.

This program reconstructs that provenance, validates every candidate, and
selects at most one candidate per foreground.  Exact selected-background loci
are unique.  Partial background/background overlaps are treated as redundancy,
not leakage: they contribute to a deterministic conflict score but are not a
hard rejection.  A positive is retained only when its paired background is
retained, keeping every emitted dataset class-balanced.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import itertools
import json
import random
import re
import shutil
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, TextIO


CHROMOSOMES = tuple(str(i) for i in range(1, 23)) + ("X",)
TRAIN_CHROMOSOMES = tuple(c for c in CHROMOSOMES if c not in ("1", "2"))
HEADER_RE = re.compile(
    r"^chr(?P<chrom>\S+)\s+(?P<start>\d+)\s+(?P<end>\d+)\s+"
    r"peak(?P<peak>\d+)_(?P<rank>\d+)_Mo_(?P<mismatch>[0-9.eE+-]+)$"
)
BIN_SIZE = 10_000


@dataclass(frozen=True)
class BedRecord:
    chrom: str
    start: int
    end: int
    text: str


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def text(self) -> str:
        return f">{self.header}\n{self.sequence}\n"


@dataclass
class Candidate:
    serial: int
    peak: int
    rank: int
    mismatch: float
    bed: BedRecord
    fasta: FastaRecord
    valid: bool = True
    rejection_reason: str = ""
    global_overlap_count: int = 0
    global_overlap_bp: int = 0


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--chromosomes",
        nargs="+",
        default=list(CHROMOSOMES),
        help="Chromosomes without the chr prefix (default: chr1--chr22 and chrX)",
    )
    parser.add_argument(
        "--output-base",
        type=Path,
        default=Path("data/processed"),
        help="Base containing task outputs and consolidation (default: data/processed)",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def bed_records(path: Path) -> Iterator[BedRecord]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty BED: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(f"{path}:{line_number}: blank BED line")
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: fewer than three columns")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid coordinates") from error
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")
            yield BedRecord(fields[0], start, end, line if line.endswith("\n") else line + "\n")


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
            else:
                if header is None:
                    if line.strip():
                        raise ValueError(f"{path}:{line_number}: sequence before header")
                    continue
                sequence.append(line.strip())
    if header is not None:
        if not sequence:
            raise ValueError(f"{path}: empty sequence for {header}")
        yield FastaRecord(header, "".join(sequence).upper())


def split_for(chromosome: str) -> str:
    return "validation" if chromosome == "1" else "test" if chromosome == "2" else "train"


def load_screen_mask(path: Path, wanted: set[str]) -> dict[str, tuple[list[int], list[tuple[int, int]]]]:
    raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for record in bed_records(path):
        if record.chrom in wanted:
            raw[record.chrom].append((record.start, record.end))
    result: dict[str, tuple[list[int], list[tuple[int, int]]]] = {}
    for chrom in wanted:
        intervals = sorted(raw.get(chrom, []))
        if not intervals:
            raise ValueError(f"SCREEN mask contains no records for {chrom}")
        merged: list[tuple[int, int]] = []
        for start, end in intervals:
            if merged and start < merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append((start, end))
        result[chrom] = ([start for start, _ in merged], merged)
    return result


def overlaps_mask(record: BedRecord, mask: tuple[list[int], list[tuple[int, int]]]) -> bool:
    starts, intervals = mask
    index = bisect.bisect_left(starts, record.end) - 1
    return index >= 0 and intervals[index][1] > record.start


def interval_overlap(first: BedRecord, second: BedRecord) -> int:
    return max(0, min(first.end, second.end) - max(first.start, second.start))


def add_global_conflict_scores(candidates: list[Candidate]) -> None:
    active: list[Candidate] = []
    for current in sorted((c for c in candidates if c.valid), key=lambda c: (c.bed.start, c.bed.end, c.serial)):
        active = [other for other in active if other.bed.end > current.bed.start]
        for other in active:
            if other.peak == current.peak:
                continue
            overlap = interval_overlap(other.bed, current.bed)
            if overlap:
                other.global_overlap_count += 1
                current.global_overlap_count += 1
                other.global_overlap_bp += overlap
                current.global_overlap_bp += overlap
        active.append(current)


def selected_conflicts(candidate: Candidate, selected: list[Candidate], bins: dict[int, set[int]]) -> tuple[int, int]:
    serials: set[int] = set()
    first_bin = candidate.bed.start // BIN_SIZE
    last_bin = (candidate.bed.end - 1) // BIN_SIZE
    for bin_number in range(first_bin, last_bin + 1):
        serials.update(bins.get(bin_number, ()))
    count = total_bp = 0
    for serial in serials:
        overlap = interval_overlap(candidate.bed, selected[serial].bed)
        if overlap:
            count += 1
            total_bp += overlap
    return count, total_bp


def add_selected(candidate: Candidate, selected: list[Candidate], bins: dict[int, set[int]]) -> None:
    selected.append(candidate)
    selected_serial = len(selected) - 1
    first_bin = candidate.bed.start // BIN_SIZE
    last_bin = (candidate.bed.end - 1) // BIN_SIZE
    for bin_number in range(first_bin, last_bin + 1):
        bins[bin_number].add(selected_serial)


def parse_candidates(
    chromosome: str,
    directory: Path,
    foregrounds: list[BedRecord],
    mask: tuple[list[int], list[tuple[int, int]]],
) -> tuple[list[Candidate], dict[int, list[Candidate]]]:
    bed_path = directory / "antinoise_candidates_2x.bed"
    fasta_path = directory / "antinoise_candidates_2x.fa"
    sentinel = object()
    candidates: list[Candidate] = []
    by_peak: dict[int, list[Candidate]] = defaultdict(list)
    streams = itertools.zip_longest(bed_records(bed_path), fasta_records(fasta_path), fillvalue=sentinel)
    for serial, pair in enumerate(streams):
        bed, fasta = pair
        if bed is sentinel or fasta is sentinel:
            raise ValueError(f"Candidate BED/FASTA mismatch for chr{chromosome}")
        assert isinstance(bed, BedRecord) and isinstance(fasta, FastaRecord)
        match = HEADER_RE.match(fasta.header)
        if match is None:
            raise ValueError(f"Malformed antinoise header: {fasta.header!r}")
        candidate = Candidate(
            serial=serial,
            peak=int(match.group("peak")),
            rank=int(match.group("rank")),
            mismatch=float(match.group("mismatch")),
            bed=bed,
            fasta=fasta,
        )
        header_coords = (f"chr{match.group('chrom')}", int(match.group("start")), int(match.group("end")))
        if header_coords != (bed.chrom, bed.start, bed.end):
            candidate.valid = False
            candidate.rejection_reason = "candidate_bed_fasta_coordinate_mismatch"
        elif candidate.peak < 0 or candidate.peak >= len(foregrounds):
            candidate.valid = False
            candidate.rejection_reason = "foreground_index_out_of_range"
        elif bed.chrom != f"chr{chromosome}":
            candidate.valid = False
            candidate.rejection_reason = "wrong_chromosome"
        elif len(fasta.sequence) != bed.end - bed.start:
            candidate.valid = False
            candidate.rejection_reason = "candidate_bed_fasta_length_mismatch"
        elif bed.end - bed.start != foregrounds[candidate.peak].end - foregrounds[candidate.peak].start:
            candidate.valid = False
            candidate.rejection_reason = "foreground_candidate_length_mismatch"
        elif set(fasta.sequence) - set("ACGT"):
            candidate.valid = False
            candidate.rejection_reason = "non_acgt_candidate_sequence"
        elif overlaps_mask(bed, mask):
            candidate.valid = False
            candidate.rejection_reason = "screen_ccre_overlap"
        candidates.append(candidate)
        by_peak[candidate.peak].append(candidate)
    for peak, items in by_peak.items():
        ranks = [candidate.rank for candidate in items]
        if len(ranks) != len(set(ranks)):
            for candidate in items:
                candidate.valid = False
                candidate.rejection_reason = "duplicate_candidate_rank"
    return candidates, by_peak


def choose_candidates(
    chromosome: str,
    foreground_count: int,
    by_peak: dict[int, list[Candidate]],
    seed: int,
) -> dict[int, Candidate]:
    eligible = {peak: [c for c in by_peak.get(peak, ()) if c.valid] for peak in range(foreground_count)}
    add_global_conflict_scores([candidate for items in eligible.values() for candidate in items])

    # One-option foregrounds are most constrained. Seeded tie-breaking avoids
    # systematically favoring the EVN-first consolidation order.
    tie_order = list(range(foreground_count))
    random.Random(seed).shuffle(tie_order)
    tie_rank = {peak: order for order, peak in enumerate(tie_order)}
    processing_order = sorted(
        range(foreground_count),
        key=lambda peak: (
            len(eligible[peak]) if eligible[peak] else 99,
            min((c.global_overlap_count for c in eligible[peak]), default=10**12),
            tie_rank[peak],
        ),
    )

    chosen: dict[int, Candidate] = {}
    exact_loci: set[tuple[str, int, int]] = set()
    selected: list[Candidate] = []
    bins: dict[int, set[int]] = defaultdict(set)
    for peak in processing_order:
        options = []
        for candidate in eligible[peak]:
            locus = (candidate.bed.chrom, candidate.bed.start, candidate.bed.end)
            if locus in exact_loci:
                continue
            overlap_count, overlap_bp = selected_conflicts(candidate, selected, bins)
            score = (
                overlap_count,
                overlap_bp,
                candidate.global_overlap_count,
                candidate.global_overlap_bp,
                candidate.rank,
                candidate.serial,
            )
            options.append((score, candidate))
        if not options:
            continue
        candidate = min(options, key=lambda item: item[0])[1]
        chosen[peak] = candidate
        exact_loci.add((candidate.bed.chrom, candidate.bed.start, candidate.bed.end))
        add_selected(candidate, selected, bins)
    return chosen


def foreground_inputs(directory: Path) -> tuple[list[BedRecord], list[FastaRecord], int]:
    evn = list(bed_records(directory / "evn_pos.bed"))
    svn = list(bed_records(directory / "svn_pos.bed"))
    consolidated = list(bed_records(directory / "consolidation.bed"))
    fasta = list(fasta_records(directory / "consolidation.fa"))
    if consolidated != evn + svn:
        raise ValueError(f"{directory}: consolidation.bed is not EVN followed by SVN")
    if len(fasta) != len(consolidated):
        raise ValueError(f"{directory}: foreground BED/FASTA count mismatch")
    for index, (bed, sequence) in enumerate(zip(consolidated, fasta)):
        if len(sequence.sequence) != bed.end - bed.start:
            raise ValueError(f"{directory}: foreground length mismatch at zero-based peak {index}")
    return consolidated, fasta, len(evn)


def write_pair(handles: dict[str, TextIO], prefix: str, foreground: BedRecord, foreground_fasta: FastaRecord, candidate: Candidate) -> None:
    handles[f"{prefix}_pos_bed"].write(foreground.text)
    handles[f"{prefix}_pos_fa"].write(foreground_fasta.text)
    handles[f"{prefix}_neg_bed"].write(candidate.bed.text)
    handles[f"{prefix}_neg_fa"].write(candidate.fasta.text)


def open_handle(path: Path) -> TextIO:
    path.parent.mkdir(parents=True, exist_ok=True)
    return path.open("w", encoding="utf-8", newline="")


def main() -> None:
    args = parse_args()
    repo = args.repo.resolve()
    chromosomes = tuple(args.chromosomes)
    unknown = sorted(set(chromosomes) - set(CHROMOSOMES))
    if unknown or len(chromosomes) != len(set(chromosomes)):
        raise ValueError(f"Invalid or repeated chromosomes: {unknown or chromosomes}")
    output_base = args.output_base if args.output_base.is_absolute() else repo / args.output_base
    consolidation_root = repo / "data/processed/consolidation/GRCh38"
    mask_path = repo / "data/raw/GRCh38-cCREs.bed"
    masks = load_screen_mask(mask_path, {f"chr{chrom}" for chrom in chromosomes})

    task_roots = {task: output_base / task / "GRCh38/real" for task in ("evn", "svn")}
    audit_root = output_base / "consolidation/GRCh38/antinoise_selection"
    destinations = [*task_roots.values(), audit_root]
    existing = [path for path in destinations if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Output exists: {existing[0]}; pass --overwrite")
    staging_parent = Path(tempfile.mkdtemp(prefix="antinoise-selection.", dir=output_base))
    staged_task_roots = {task: staging_parent / task / "GRCh38/real" for task in ("evn", "svn")}
    staged_audit_root = staging_parent / "consolidation/GRCh38/antinoise_selection"

    handles: dict[str, TextIO] = {}
    counts: dict[str, dict[str, int]] = {task: defaultdict(int) for task in ("evn", "svn")}
    summary_rows: list[dict[str, object]] = []
    try:
        for task, root in staged_task_roots.items():
            for split in ("validation", "test", "train_100"):
                for cls in ("pos", "neg"):
                    for extension in ("bed", "fa"):
                        handles[f"{task}_{split}_{cls}_{extension}"] = open_handle(root / f"{split}_{cls}.{extension}")
        staged_audit_root.mkdir(parents=True, exist_ok=True)
        provenance_handle = open_handle(staged_audit_root / "provenance.tsv")
        provenance = csv.writer(provenance_handle, delimiter="\t", lineterminator="\n")
        provenance.writerow((
            "chromosome", "split", "task", "peak_id_zero_based", "foreground_row_one_based",
            "foreground_start", "foreground_end", "foreground_id", "candidate_count",
            "valid_candidate_count", "selected_rank", "negative_start", "negative_end",
            "at_mismatch", "global_candidate_overlap_count", "global_candidate_overlap_bp", "status",
        ))

        for chromosome in chromosomes:
            directory = consolidation_root / "chromosomes" / f"chrom{chromosome}"
            foregrounds, foreground_fastas, evn_boundary = foreground_inputs(directory)
            candidates, by_peak = parse_candidates(
                chromosome, directory, foregrounds, masks[f"chr{chromosome}"]
            )
            chosen = choose_candidates(chromosome, len(foregrounds), by_peak, args.seed + CHROMOSOMES.index(chromosome))
            split = split_for(chromosome)
            chromosome_counts = {"evn": 0, "svn": 0}
            for peak, foreground in enumerate(foregrounds):
                task = "evn" if peak < evn_boundary else "svn"
                candidate = chosen.get(peak)
                fields = foreground.text.rstrip("\n").split("\t")
                foreground_id = fields[3] if len(fields) > 3 else ""
                valid_count = sum(candidate_item.valid for candidate_item in by_peak.get(peak, ()))
                if candidate is None:
                    provenance.writerow((
                        f"chr{chromosome}", split, task, peak, peak + 1, foreground.start,
                        foreground.end, foreground_id, len(by_peak.get(peak, ())), valid_count,
                        "", "", "", "", "", "", "no_acceptable_candidate",
                    ))
                    continue
                destination = "train_100" if split == "train" else split
                write_pair(
                    handles,
                    f"{task}_{destination}",
                    foreground,
                    foreground_fastas[peak],
                    candidate,
                )
                if split == "train":
                    chromosome_root = staged_task_roots[task] / "chromosomes" / f"chrom{chromosome}"
                    chromosome_handles = {}
                    for cls in ("pos", "neg"):
                        for extension in ("bed", "fa"):
                            key = f"{task}_chrom{chromosome}_{cls}_{extension}"
                            if key not in handles:
                                handles[key] = open_handle(chromosome_root / f"train_{cls}.{extension}")
                            chromosome_handles[f"chrom_{cls}_{extension}"] = handles[key]
                    write_pair(
                        chromosome_handles,
                        "chrom",
                        foreground,
                        foreground_fastas[peak],
                        candidate,
                    )
                chromosome_counts[task] += 1
                counts[task][split] += 1
                provenance.writerow((
                    f"chr{chromosome}", split, task, peak, peak + 1, foreground.start,
                    foreground.end, foreground_id, len(by_peak.get(peak, ())), valid_count,
                    candidate.rank, candidate.bed.start, candidate.bed.end, candidate.mismatch,
                    candidate.global_overlap_count, candidate.global_overlap_bp, "selected",
                ))
            summary_rows.append({
                "chromosome": f"chr{chromosome}",
                "split": split,
                "foregrounds": len(foregrounds),
                "candidates": len(candidates),
                "selected_evn": chromosome_counts["evn"],
                "selected_svn": chromosome_counts["svn"],
                "unresolved": len(foregrounds) - sum(chromosome_counts.values()),
            })

        for handle in handles.values():
            handle.close()
        provenance_handle.close()
        with (staged_audit_root / "summary.json").open("w", encoding="utf-8") as handle:
            json.dump(
                {
                    "seed": args.seed,
                    "chromosomes": [f"chr{chrom}" for chrom in chromosomes],
                    "peak_indexing": "zero-based; BED/FASTA foreground row = peak_id + 1",
                    "positive_negative_overlap_rule": "reject any overlap with SCREEN cCREs",
                    "negative_duplicate_rule": "reject exact selected-locus reuse within chromosome",
                    "partial_negative_overlap_rule": "allowed; minimized by deterministic greedy score",
                    "records_per_class": {task: dict(values) for task, values in counts.items()},
                    "chromosome_results": summary_rows,
                },
                handle,
                indent=2,
                sort_keys=True,
            )
            handle.write("\n")

        for destination in destinations:
            relative = destination.relative_to(output_base)
            source = staging_parent / relative
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
        shutil.rmtree(staging_parent, ignore_errors=True)

    for task in ("evn", "svn"):
        values = counts[task]
        print(
            f"{task.upper()}: validation={values['validation']:,}; "
            f"test={values['test']:,}; training={values['train']:,}"
        )
    print(f"Audit: {audit_root}")


if __name__ == "__main__":
    main()
