#!/usr/bin/env python3
"""Download JASPAR 2026 CORE motifs and build a human-relevant MEME file.

The JASPAR CORE non-redundant MEME file contains motifs from several taxa.
This script retains the current non-redundant motif for every JASPAR base-ID
lineage that has at least one Homo sapiens (NCBI taxonomy ID 9606) metadata
record.  Matching by base ID intentionally carries human evidence between
versions of the same JASPAR motif lineage.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Iterator


MEME_URL = (
    "https://jaspar.elixir.no/download/data/2026/CORE/"
    "JASPAR2026_CORE_non-redundant_pfms_meme.txt"
)
METADATA_URL = (
    "https://mencius-pub.uio.no/JASPAR/JASPAR_metadata/2026/"
    "ultimate_metadata_table_CORE.tsv"
)
MEME_NAME = "JASPAR2026_CORE_non-redundant_pfms_meme.txt"
METADATA_NAME = "ultimate_metadata_table_CORE.tsv"
HUMAN_MEME_NAME = "JASPAR2026_CORE_non-redundant_pfms_homo_sapiens_meme.txt"
AUDIT_NAME = "JASPAR2026_CORE_non-redundant_pfms_homo_sapiens_audit.tsv"
SUMMARY_NAME = "JASPAR2026_CORE_non-redundant_pfms_homo_sapiens_summary.json"
CHECKSUM_NAME = "checksums.sha256"
HUMAN_TAX_ID = "9606"


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "data/motifs",
        help="Destination for downloaded sources and derived outputs",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "ccre-classification/JASPAR2026-builder"},
    )
    with urllib.request.urlopen(request) as response, destination.open("wb") as output:
        while chunk := response.read(1024 * 1024):
            output.write(chunk)
        output.flush()
        os.fsync(output.fileno())
    if destination.stat().st_size == 0:
        raise RuntimeError(f"Downloaded an empty file from {url}")


def read_metadata(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"matrix_id", "base_id", "version", "tax_id", "species"}
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Metadata is missing columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("Metadata table has no records")
    return rows


def meme_records(path: Path) -> tuple[list[str], list[tuple[str, str, list[str]]]]:
    preamble: list[str] = []
    records: list[tuple[str, str, list[str]]] = []
    block: list[str] = []
    matrix_id = ""
    motif_name = ""
    with path.open(encoding="utf-8", newline="") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("MOTIF "):
                if block:
                    records.append((matrix_id, motif_name, block))
                fields = line.rstrip("\r\n").split(maxsplit=2)
                if len(fields) < 2:
                    raise ValueError(f"{path}:{line_number}: malformed MOTIF header")
                matrix_id = fields[1]
                motif_name = fields[2] if len(fields) == 3 else ""
                block = [line]
            elif block:
                block.append(line)
            else:
                preamble.append(line)
    if block:
        records.append((matrix_id, motif_name, block))
    if not records:
        raise ValueError(f"No MOTIF records found in {path}")
    ids = [record[0] for record in records]
    if len(ids) != len(set(ids)):
        raise ValueError("MEME file contains duplicate matrix IDs")
    return preamble, records


def base_id(matrix_id: str) -> str:
    if "." not in matrix_id:
        raise ValueError(f"JASPAR matrix ID has no version suffix: {matrix_id}")
    base, version = matrix_id.rsplit(".", 1)
    if not version.isdigit():
        raise ValueError(f"Invalid JASPAR matrix ID: {matrix_id}")
    return base


def write_filtered_meme(
    destination: Path,
    preamble: list[str],
    records: list[tuple[str, str, list[str]]],
    human_base_ids: set[str],
) -> int:
    retained = 0
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(preamble)
        for matrix_id, _, block in records:
            if base_id(matrix_id) in human_base_ids:
                handle.writelines(block)
                retained += 1
        handle.flush()
        os.fsync(handle.fileno())
    return retained


def atomic_targets(output_dir: Path) -> Iterator[tuple[Path, Path]]:
    """Yield final and temporary paths, retaining temporaries until promotion."""
    names = (MEME_NAME, METADATA_NAME, HUMAN_MEME_NAME, AUDIT_NAME, SUMMARY_NAME, CHECKSUM_NAME)
    for name in names:
        handle = tempfile.NamedTemporaryFile(
            dir=output_dir, prefix=f".{name}.", suffix=".tmp", delete=False
        )
        temporary = Path(handle.name)
        handle.close()
        yield output_dir / name, temporary


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs = list(atomic_targets(output_dir))
    paths = {final.name: temporary for final, temporary in pairs}
    existing = [final for final, _ in pairs if final.exists()]
    if existing and not args.overwrite:
        for _, temporary in pairs:
            temporary.unlink(missing_ok=True)
        raise FileExistsError(f"Output exists: {existing[0]}; pass --overwrite")

    try:
        download(MEME_URL, paths[MEME_NAME])
        download(METADATA_URL, paths[METADATA_NAME])
        metadata = read_metadata(paths[METADATA_NAME])
        preamble, records = meme_records(paths[MEME_NAME])

        metadata_matrix_ids = {row["matrix_id"] for row in metadata}
        meme_matrix_ids = {record[0] for record in records}
        missing_metadata = meme_matrix_ids.difference(metadata_matrix_ids)
        if missing_metadata:
            preview = ", ".join(sorted(missing_metadata)[:5])
            raise ValueError(f"MEME IDs absent from metadata: {preview}")

        human_rows = [row for row in metadata if row["tax_id"] == HUMAN_TAX_ID]
        if not human_rows:
            raise ValueError("No Homo sapiens metadata rows (tax_id 9606)")
        if any(row["species"] != "Homo sapiens" for row in human_rows):
            raise ValueError("tax_id 9606 is not consistently labelled Homo sapiens")
        human_matrix_ids = {row["matrix_id"] for row in human_rows}
        human_by_base: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in human_rows:
            human_by_base[row["base_id"]].append(row)
        human_base_ids = set(human_by_base)

        retained = write_filtered_meme(
            paths[HUMAN_MEME_NAME], preamble, records, human_base_ids
        )
        exact_matches = 0
        lineage_matches = 0
        with paths[AUDIT_NAME].open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                (
                    "matrix_id",
                    "base_id",
                    "motif_name",
                    "retained",
                    "human_match_type",
                    "human_evidence_matrix_ids",
                )
            )
            for matrix_id, motif_name, _ in records:
                base = base_id(matrix_id)
                keep = base in human_base_ids
                if matrix_id in human_matrix_ids:
                    match_type = "exact_matrix_id"
                    exact_matches += 1
                elif keep:
                    match_type = "base_id_lineage"
                    lineage_matches += 1
                else:
                    match_type = "none"
                evidence = ";".join(
                    sorted(row["matrix_id"] for row in human_by_base.get(base, ()))
                )
                writer.writerow((matrix_id, base, motif_name, str(keep).lower(), match_type, evidence))
            handle.flush()
            os.fsync(handle.fileno())

        if retained != exact_matches + lineage_matches:
            raise AssertionError("Retained motif accounting mismatch")
        summary = {
            "release": "JASPAR 2026",
            "collection": "CORE",
            "source_meme_url": MEME_URL,
            "source_metadata_url": METADATA_URL,
            "filter": {
                "tax_id": HUMAN_TAX_ID,
                "species": "Homo sapiens",
                "matching_key": "base_id",
                "rationale": (
                    "Retain the current non-redundant motif version when any version "
                    "in its JASPAR base-ID lineage has Homo sapiens evidence."
                ),
            },
            "counts": {
                "source_motifs": len(records),
                "source_metadata_rows": len(metadata),
                "human_metadata_rows": len(human_rows),
                "retained_exact_matrix_id": exact_matches,
                "retained_base_id_lineage": lineage_matches,
                "retained_total": retained,
            },
        }
        paths[SUMMARY_NAME].write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        checksummed = (MEME_NAME, METADATA_NAME, HUMAN_MEME_NAME, AUDIT_NAME, SUMMARY_NAME)
        with paths[CHECKSUM_NAME].open("w", encoding="utf-8", newline="") as handle:
            for name in checksummed:
                handle.write(f"{sha256(paths[name])}  {name}\n")
            handle.flush()
            os.fsync(handle.fileno())

        for final, temporary in pairs:
            os.chmod(temporary, 0o644)
            os.replace(temporary, final)
    except BaseException:
        for _, temporary in pairs:
            temporary.unlink(missing_ok=True)
        raise

    print(f"Source motifs: {len(records):,}")
    print(f"Exact Homo sapiens versions retained: {exact_matches:,}")
    print(f"Base-ID lineage versions retained: {lineage_matches:,}")
    print(f"Human-relevant motifs written: {retained:,}")
    print(f"Output: {output_dir / HUMAN_MEME_NAME}")


if __name__ == "__main__":
    main()
