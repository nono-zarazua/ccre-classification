#!/usr/bin/env python3
"""Build a MEME database from the HOCOMOCO v14 clustered CORE motifs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path


BASE_URL = (
    "https://hocomoco14.autosome.org/final_bundle/hocomoco14/"
    "H14CORE-CLUSTERED"
)
PFM_NAME = "H14CORE-CLUSTERED_pfm.tar.gz"
CLUSTERS_NAME = "cluster_list.tsv"
ANNOTATION_NAME = "H14CORE-CLUSTERED_annotation.jsonl"
MEME_NAME = "H14CORE-CLUSTERED.meme"
SUMMARY_NAME = "H14CORE-CLUSTERED_summary.json"
CHECKSUM_NAME = "H14CORE-CLUSTERED_checksums.sha256"


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=repo / "data/motifs/hocomoco14")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def download(url: str, destination: Path) -> None:
    request = urllib.request.Request(
        url, headers={"User-Agent": "ccre-classification/HOCOMOCO14-builder"}
    )
    with urllib.request.urlopen(request) as response, destination.open("wb") as out:
        shutil.copyfileobj(response, out)
        out.flush()
        os.fsync(out.fileno())
    if destination.stat().st_size == 0:
        raise RuntimeError(f"Downloaded an empty file from {url}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_clusters(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "Representative_Motif", "Cluster_Size", "Clustered_Motifs",
            "Primary_Family", "Primary_Subfamily",
            "Representative_Motif_Gene_Symbol", "Clustered_Motifs_Gene_Symbols",
        }
        missing = required.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"Cluster table is missing columns: {sorted(missing)}")
        rows = list(reader)
    if len(rows) != 648:
        raise ValueError(f"Expected 648 H14CORE clusters, found {len(rows)}")
    representatives = [row["Representative_Motif"] for row in rows]
    if len(representatives) != len(set(representatives)):
        raise ValueError("Cluster table contains duplicate representative motifs")
    return rows


def read_annotations(path: Path) -> dict[str, dict]:
    records: dict[str, dict] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            record = json.loads(line)
            name = record.get("name")
            if not name or name in records:
                raise ValueError(f"Invalid or duplicate annotation at line {line_number}")
            records[name] = record
    if len(records) != 648:
        raise ValueError(f"Expected 648 representative annotations, found {len(records)}")
    return records


def read_pfms(archive_path: Path) -> dict[str, list[list[float]]]:
    matrices: dict[str, list[list[float]]] = {}
    with tarfile.open(archive_path, "r:gz") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        for member in members:
            extracted = archive.extractfile(member)
            if extracted is None:
                raise RuntimeError(f"Cannot read {member.name}")
            lines = extracted.read().decode("utf-8").strip().splitlines()
            if not lines or not lines[0].startswith(">"):
                raise ValueError(f"Malformed PFM: {member.name}")
            name = lines[0][1:].strip()
            matrix = [[float(value) for value in line.split()] for line in lines[1:]]
            if not matrix or any(len(row) != 4 for row in matrix):
                raise ValueError(f"PFM must contain L rows of A/C/G/T probabilities: {name}")
            for row in matrix:
                if any(value < 0 for value in row) or not 0.999 <= sum(row) <= 1.001:
                    raise ValueError(f"Invalid probability row in {name}: {row}")
            if name in matrices:
                raise ValueError(f"Duplicate PFM: {name}")
            matrices[name] = matrix
    if len(matrices) != 648:
        raise ValueError(f"Expected 648 representative PFMs, found {len(matrices)}")
    return matrices


def write_meme(path: Path, clusters: list[dict[str, str]], matrices: dict[str, list[list[float]]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(
            "MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
            "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n"
        )
        for cluster in clusters:
            representative = cluster["Representative_Motif"]
            matrix = matrices[representative]
            handle.write(f"MOTIF {representative}\n")
            handle.write(
                f"letter-probability matrix: alength= 4 w= {len(matrix)} "
                "nsites= 1 E= 0\n"
            )
            for row in matrix:
                handle.write(" ".join(f"{value:.16g}" for value in row) + "\n")
            handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    names = (PFM_NAME, CLUSTERS_NAME, ANNOTATION_NAME, MEME_NAME, SUMMARY_NAME, CHECKSUM_NAME)
    existing = [output_dir / name for name in names if (output_dir / name).exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Output exists: {existing[0]}; pass --overwrite")

    work = Path(tempfile.mkdtemp(prefix=".hocomoco14.", dir=output_dir))
    try:
        pfm = work / PFM_NAME
        clusters_path = work / CLUSTERS_NAME
        annotation_path = work / ANNOTATION_NAME
        download(f"{BASE_URL}/{PFM_NAME}", pfm)
        download(f"{BASE_URL}/{CLUSTERS_NAME}", clusters_path)
        download(f"{BASE_URL}/{ANNOTATION_NAME}", annotation_path)
        clusters = read_clusters(clusters_path)
        annotations = read_annotations(annotation_path)
        matrices = read_pfms(pfm)
        representatives = {row["Representative_Motif"] for row in clusters}
        if representatives != matrices.keys() or representatives != annotations.keys():
            raise ValueError("PFM, cluster-table, and annotation identifiers do not agree")

        meme = work / MEME_NAME
        write_meme(meme, clusters, matrices)
        summary = {
            "release": "HOCOMOCO v14",
            "collection": "H14CORE-CLUSTERED",
            "motif_clusters": len(clusters),
            "alphabet_order": "ACGT",
            "source_urls": {
                PFM_NAME: f"{BASE_URL}/{PFM_NAME}",
                CLUSTERS_NAME: f"{BASE_URL}/{CLUSTERS_NAME}",
                ANNOTATION_NAME: f"{BASE_URL}/{ANNOTATION_NAME}",
            },
            "purpose": (
                "Clustered representative motifs prevent near-duplicate TF PWMs "
                "from fragmenting motif counts. Cluster-level labels do not identify "
                "one uniquely active transcription factor."
            ),
        }
        (work / SUMMARY_NAME).write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        with (work / CHECKSUM_NAME).open("w", encoding="utf-8") as handle:
            for name in names[:-1]:
                handle.write(f"{sha256(work / name)}  {name}\n")

        for name in names:
            os.replace(work / name, output_dir / name)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
