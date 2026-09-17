"""Shared BED/FASTA validation and sampling helpers for GRCh38 training sets."""

from __future__ import annotations

import itertools
import random
import sqlite3
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO




@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    text: str


@dataclass(frozen=True)
class SourcePair:
    bed: Path
    fasta: Path
    count: int



def fasta_records(path: Path) -> Iterator[FastaRecord]:
    header: str | None = None
    lines: list[str] = []
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(">"):
                if header is not None:
                    if len(lines) == 1:
                        raise ValueError(f"{path}: {header!r} has no sequence")
                    yield FastaRecord(header, "".join(lines))
                header = line[1:].strip()
                if not header:
                    raise ValueError(f"{path}:{line_number}: empty FASTA header")
                lines = [line]
            elif header is None:
                if line.strip():
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
            else:
                lines.append(line)
    if header is not None:
        if len(lines) == 1:
            raise ValueError(f"{path}: {header!r} has no sequence")
        text = "".join(lines)
        yield FastaRecord(header, text if text.endswith("\n") else text + "\n")


def nonempty_bed_lines(path: Path) -> Iterator[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                raise ValueError(f"{path}:{line_number}: blank BED record")
            yield line if line.endswith("\n") else line + "\n"


def count_bed(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty BED: {path}")
    return sum(1 for _ in nonempty_bed_lines(path))


def count_fasta(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty FASTA: {path}")
    return sum(1 for _ in fasta_records(path))




def unique_coordinate_candidates(
    source: SourcePair, already_selected: set[tuple[str, int, int]]
) -> list[tuple[int, tuple[str, int, int]]]:
    """Return one candidate per exact locus, excluding loci selected previously."""
    candidates: list[tuple[int, tuple[str, int, int]]] = []
    local_seen: set[tuple[str, int, int]] = set()
    for index, line in enumerate(nonempty_bed_lines(source.bed)):
        fields = line.split()
        if len(fields) < 3:
            raise ValueError(f"Malformed BED record in {source.bed}: {line.rstrip()}")
        key = (fields[0], int(fields[1]), int(fields[2]))
        if key in local_seen or key in already_selected:
            continue
        local_seen.add(key)
        candidates.append((index, key))
    return candidates


def sample_candidates(
    candidates: list[tuple[int, tuple[str, int, int]]], quota: int, seed: int
) -> tuple[set[int], set[tuple[str, int, int]]]:
    if quota > len(candidates):
        raise ValueError(f"Requested {quota:,} records from {len(candidates):,} unique loci")
    chosen = candidates if quota == len(candidates) else random.Random(seed).sample(candidates, quota)
    return {index for index, _ in chosen}, {key for _, key in chosen}


def copy_selected_pair(
    source: SourcePair,
    indices: set[int],
    bed_output: TextIO,
    fasta_output: TextIO,
) -> int:
    written = 0
    sentinel = object()
    pairs = itertools.zip_longest(
        nonempty_bed_lines(source.bed), fasta_records(source.fasta), fillvalue=sentinel
    )
    for index, pair in enumerate(pairs):
        bed_record, fasta_record = pair
        if bed_record is sentinel or fasta_record is sentinel:
            raise ValueError(f"BED/FASTA stream mismatch: {source.bed} and {source.fasta}")
        if index in indices:
            bed_output.write(bed_record)
            fasta_output.write(fasta_record.text)
            written += 1
    if written != len(indices):
        raise AssertionError(f"Expected {len(indices):,} selected records; wrote {written:,}")
    return written


def validate_unique_fasta(path: Path) -> int:
    database = tempfile.NamedTemporaryFile(prefix="fasta_ids.", suffix=".sqlite", delete=False)
    database_path = Path(database.name)
    database.close()
    try:
        connection = sqlite3.connect(database_path)
        connection.execute("CREATE TABLE ids (id TEXT PRIMARY KEY)")
        count = 0
        for record in fasta_records(path):
            try:
                connection.execute("INSERT INTO ids VALUES (?)", (record.identifier,))
            except sqlite3.IntegrityError as error:
                raise ValueError(f"Duplicate FASTA identifier in {path}: {record.identifier}") from error
            count += 1
        connection.close()
        return count
    finally:
        database_path.unlink(missing_ok=True)


def validate_unique_bed(paths: list[Path], label: str) -> int:
    """Verify exact genomic loci are unique within and across all supplied splits."""
    seen: set[tuple[str, int, int]] = set()
    count = 0
    for path in paths:
        for line in nonempty_bed_lines(path):
            fields = line.split()
            key = (fields[0], int(fields[1]), int(fields[2]))
            if key in seen:
                raise ValueError(f"Duplicate exact locus across {label}: {key}")
            seen.add(key)
            count += 1
    return count


