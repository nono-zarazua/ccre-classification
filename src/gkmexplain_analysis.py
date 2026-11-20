#!/usr/bin/env python3
"""Prepare null sequences and summarize final-model GkmExplain outputs."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import random
import shutil
import sys
import tempfile
import time
from collections import defaultdict
from itertools import zip_longest
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from deeplift.dinuc_shuffle import dinuc_shuffle
from memelite import tomtom
from memelite.io import read_meme
from modiscolite import io as modisco_io
from modiscolite import tfmodisco


DNA = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    null_parser = subparsers.add_parser("make-null")
    null_parser.add_argument("--input-fasta", type=Path, required=True)
    null_parser.add_argument("--output-fasta", type=Path, required=True)
    null_parser.add_argument("--n-sequences", type=int, default=1000)
    null_parser.add_argument("--seed", type=int, default=42)

    summary = subparsers.add_parser("summarize")
    summary.add_argument("--task", choices=("evn", "svn", "evs"), required=True)
    summary.add_argument("--model", type=Path, required=True)
    summary.add_argument("--positive-fasta", type=Path, required=True)
    summary.add_argument("--negative-fasta", type=Path, required=True)
    summary.add_argument("--raw-dir", type=Path, required=True)
    summary.add_argument("--motif-meme", type=Path, required=True)
    summary.add_argument("--output-dir", type=Path, required=True)
    summary.add_argument("--n-jobs", type=int, default=1)
    summary.add_argument(
        "--score-shards",
        type=int,
        default=8,
        help="Number of round-robin shards used to create GkmExplain outputs",
    )
    summary.add_argument("--kmer-length", type=int, default=0)
    summary.add_argument("--window-size", type=int, default=21)
    summary.add_argument("--target-seqlet-fdr", type=float, default=0.15)
    summary.add_argument("--tomtom-pvalue", type=float, default=0.05)
    summary.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    parts: list[str] = []
    with path.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(parts).upper()))
                identifier = line[1:].strip()
                if not identifier:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                parts = []
            else:
                if identifier is None:
                    raise ValueError(f"{path}:{line_number}: sequence before header")
                parts.append(line)
    if identifier is not None:
        records.append((identifier, "".join(parts).upper()))
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    identifiers = [identifier for identifier, _ in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate FASTA identifiers in {path}")
    for identifier, sequence in records:
        unknown = sorted(set(sequence) - set(DNA))
        if unknown:
            raise ValueError(f"Non-ACGT bases {unknown} in {path}: {identifier}")
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n{sequence}\n")


def make_null(args: argparse.Namespace) -> None:
    if args.n_sequences <= 0:
        raise ValueError("--n-sequences must be positive")
    if args.output_fasta.is_file() and args.output_fasta.stat().st_size > 0:
        print(f"SKIP: {args.output_fasta}")
        return
    records = read_fasta(args.input_fasta)
    n = min(args.n_sequences, len(records))
    rng = np.random.RandomState(args.seed)
    selected = np.sort(rng.choice(len(records), size=n, replace=False))
    shuffled = []
    for index in selected:
        identifier, sequence = records[index]
        shuffled_one_hot = dinuc_shuffle(one_hot(sequence), rng=rng)
        shuffled_sequence = "".join(DNA[int(row.argmax())] for row in shuffled_one_hot)
        shuffled.append((f"dinuc_shuffle::{identifier}", shuffled_sequence))
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output_fasta.with_name(f".{args.output_fasta.name}.tmp")
    write_fasta(temporary, shuffled)
    temporary.replace(args.output_fasta)
    print(f"Wrote {n:,} dinucleotide-shuffled sequences to {args.output_fasta}")


def sharded_record_order(
    records: list[tuple[str, str]], n_shards: int
) -> list[tuple[str, str]]:
    """Reconstruct the shard-concatenated order used by the scoring wrapper."""
    if n_shards <= 0:
        raise ValueError("Number of GkmExplain score shards must be positive")
    return [
        record
        for shard in range(n_shards)
        for record in records[shard::n_shards]
    ]


def parse_explanations(
    path: Path,
    records: list[tuple[str, str]],
    n_shards: int,
) -> dict[str, tuple[float, np.ndarray]]:
    """Parse scores while restoring full FASTA headers lost by GkmExplain.

    GkmExplain reports only the first whitespace-delimited token from each FASTA
    header. The parallel scoring wrapper distributes records round-robin and
    concatenates outputs shard-by-shard, so the full identifiers can be restored
    deterministically from the original FASTA and shard count.
    """
    explanations: dict[str, tuple[float, np.ndarray]] = {}
    expected_records = sharded_record_order(records, n_shards)
    with path.open() as handle:
        for line_number, pair in enumerate(
            zip_longest(handle, expected_records), start=1
        ):
            raw, expected = pair
            if raw is None or expected is None:
                raise ValueError(
                    f"{path}: row count does not match {len(expected_records):,} "
                    "FASTA records"
                )
            fields = raw.rstrip("\n").split("\t", maxsplit=2)
            if len(fields) != 3:
                raise ValueError(f"{path}:{line_number}: expected three tab-separated fields")
            reported_identifier, prediction, encoded = fields
            identifier, sequence = expected
            expected_reported_identifier = identifier.split(maxsplit=1)[0]
            if reported_identifier != expected_reported_identifier:
                raise ValueError(
                    f"{path}:{line_number}: reported identifier "
                    f"{reported_identifier!r} does not match expected FASTA token "
                    f"{expected_reported_identifier!r}"
                )
            if identifier in explanations:
                raise ValueError(f"Duplicate full FASTA identifier in {path}: {identifier}")
            rows = [row for row in encoded.split(";") if row]
            matrix = np.asarray(
                [[float(value) for value in row.split(",")] for row in rows],
                dtype=np.float32,
            )
            if matrix.ndim != 2 or matrix.shape[1] != 4:
                raise ValueError(f"{path}:{line_number}: invalid explanation shape {matrix.shape}")
            if matrix.shape[0] != len(sequence):
                raise ValueError(
                    f"{path}:{line_number}: explanation length {matrix.shape[0]} "
                    f"does not match sequence length {len(sequence)}"
                )
            explanations[identifier] = (float(prediction), matrix)
    if not explanations:
        raise ValueError(f"No explanations in {path}")
    return explanations


def one_hot(sequence: str) -> np.ndarray:
    result = np.zeros((len(sequence), 4), dtype=np.int8)
    for index, base in enumerate(sequence):
        result[index, DNA.index(base)] = 1
    return result


def normalize_hypothetical(
    observed: np.ndarray, hypothetical: np.ndarray
) -> np.ndarray:
    observed_per_position = observed.sum(axis=1)
    sign = np.sign(observed_per_position)
    same_sign = np.sign(hypothetical) * sign[:, None] > 0
    denominator = (hypothetical * same_sign).sum(axis=1)
    ratio = np.divide(
        observed_per_position,
        denominator,
        out=np.zeros_like(observed_per_position),
        where=denominator != 0,
    )
    normalized = hypothetical * ratio[:, None]
    if not np.isfinite(normalized).all():
        raise RuntimeError("Hypothetical-score normalization produced non-finite values")
    return normalized.astype(np.float32, copy=False)


def load_class_arrays(
    records: list[tuple[str, str]],
    observed_path: Path,
    hypothetical_path: Path,
    score_shards: int,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]], list[np.ndarray]]:
    observed = parse_explanations(observed_path, records, score_shards)
    hypothetical = parse_explanations(hypothetical_path, records, score_shards)
    expected_ids = {identifier for identifier, _ in records}
    if set(observed) != expected_ids or set(hypothetical) != expected_ids:
        raise RuntimeError("FASTA and GkmExplain identifiers are not synchronized")
    max_length = max(len(sequence) for _, sequence in records)
    X = np.zeros((len(records), max_length, 4), dtype=np.int8)
    H = np.zeros((len(records), max_length, 4), dtype=np.float32)
    manifest: list[dict[str, Any]] = []
    per_position: list[np.ndarray] = []
    for example_idx, (identifier, sequence) in enumerate(records):
        prediction, observed_matrix = observed[identifier]
        hypothetical_prediction, hypothetical_matrix = hypothetical[identifier]
        if observed_matrix.shape[0] != len(sequence):
            raise RuntimeError(f"Observed explanation length mismatch: {identifier}")
        if hypothetical_matrix.shape != observed_matrix.shape:
            raise RuntimeError(f"Hypothetical explanation shape mismatch: {identifier}")
        if not np.isclose(prediction, hypothetical_prediction, rtol=1e-6, atol=1e-8):
            raise RuntimeError(f"Prediction differs between modes for {identifier}")
        encoded = one_hot(sequence)
        normalized = normalize_hypothetical(observed_matrix, hypothetical_matrix)
        left = (max_length - len(sequence)) // 2
        right = left + len(sequence)
        X[example_idx, left:right] = encoded
        H[example_idx, left:right] = normalized
        per_position.append(observed_matrix.sum(axis=1))
        manifest.append(
            {
                "example_idx": example_idx,
                "sequence_id": identifier,
                "length": len(sequence),
                "left_padding": left,
                "prediction": prediction,
            }
        )
    return X, H, manifest, per_position


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def canonical(sequence: str) -> str:
    reverse = reverse_complement(sequence)
    return min(sequence, reverse)


def accumulate_kmers(
    records: list[tuple[str, str]],
    scores: list[np.ndarray],
    class_index: int,
    length: int,
    accumulators: dict[str, list[float]],
) -> None:
    offset = class_index * 5
    for sequence_index, ((_, sequence), position_scores) in enumerate(zip(records, scores)):
        if len(sequence) < length:
            continue
        prefix = np.concatenate(([0.0], np.cumsum(position_scores, dtype=np.float64)))
        seen: set[str] = set()
        for start in range(len(sequence) - length + 1):
            kmer = canonical(sequence[start : start + length])
            value = float(prefix[start + length] - prefix[start])
            stats = accumulators[kmer]
            stats[offset] += 1
            stats[offset + 2] += value
            stats[offset + 3] += abs(value)
            if kmer not in seen:
                stats[offset + 1] += 1
                seen.add(kmer)


def kmer_table(accumulators: dict[str, list[float]]) -> pd.DataFrame:
    rows = []
    for kmer, values in accumulators.items():
        pos_occ, pos_seq, pos_sum, pos_abs, _, neg_occ, neg_seq, neg_sum, neg_abs, _ = values
        rows.append(
            {
                "canonical_kmer": kmer,
                "positive_occurrences": int(pos_occ),
                "positive_sequences": int(pos_seq),
                "positive_total_contribution": pos_sum,
                "positive_mean_contribution": pos_sum / pos_occ if pos_occ else np.nan,
                "positive_mean_absolute_contribution": pos_abs / pos_occ if pos_occ else np.nan,
                "negative_occurrences": int(neg_occ),
                "negative_sequences": int(neg_seq),
                "negative_total_contribution": neg_sum,
                "negative_mean_contribution": neg_sum / neg_occ if neg_occ else np.nan,
                "negative_mean_absolute_contribution": neg_abs / neg_occ if neg_occ else np.nan,
                "mean_contribution_contrast": (
                    (pos_sum / pos_occ if pos_occ else 0.0)
                    - (neg_sum / neg_occ if neg_occ else 0.0)
                ),
            }
        )
    table = pd.DataFrame(rows)
    return table.sort_values(
        ["positive_mean_absolute_contribution", "negative_mean_absolute_contribution"],
        ascending=False,
    ).reset_index(drop=True)


def bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    if len(pvalues) == 0:
        return np.asarray([], dtype=float)
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    adjusted = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def annotate_patterns(
    class_name: str,
    pos_patterns: Any,
    neg_patterns: Any,
    motifs: dict[str, np.ndarray],
    n_jobs: int,
    threshold: float,
) -> pd.DataFrame:
    pattern_records: list[tuple[str, int, Any]] = []
    for sign, patterns in (("positive", pos_patterns), ("negative", neg_patterns)):
        for index, pattern in enumerate(patterns or []):
            pattern_records.append((sign, index, pattern))
    if not pattern_records:
        return pd.DataFrame(
            columns=[
                "class_name", "sign", "pattern_index", "n_seqlets",
                "pattern_length", "jaspar_motif", "tomtom_pvalue",
                "tomtom_qvalue", "significant_match",
            ]
        )
    queries = [np.asarray(pattern.sequence).T for _, _, pattern in pattern_records]
    motif_names = list(motifs)
    result = tomtom(queries, list(motifs.values()), n_nearest=1, n_jobs=n_jobs)
    pvalues = result[0, :, 0]
    motif_indices = result[5, :, 0].astype(int)
    qvalues = bh_adjust(pvalues)
    rows = []
    for record_index, (sign, pattern_index, pattern) in enumerate(pattern_records):
        rows.append(
            {
                "class_name": class_name,
                "sign": sign,
                "pattern_index": pattern_index,
                "n_seqlets": len(pattern.seqlets),
                "pattern_length": int(pattern.sequence.shape[0]),
                "total_contribution": float(np.asarray(pattern.contrib_scores).sum()),
                "jaspar_motif": motif_names[motif_indices[record_index]],
                "tomtom_pvalue": float(pvalues[record_index]),
                "tomtom_qvalue": float(qvalues[record_index]),
                "significant_match": bool(qvalues[record_index] <= threshold),
            }
        )
    return pd.DataFrame(rows)


def window_diagnostics(
    class_name: str,
    source: str,
    records: list[tuple[str, str]],
    explanation_path: Path,
    window_size: int,
    score_shards: int,
) -> pd.DataFrame:
    explanations = parse_explanations(explanation_path, records, score_shards)
    rows = []
    for identifier, sequence in records:
        prediction, matrix = explanations[identifier]
        scores = matrix.sum(axis=1)
        width = min(window_size, len(scores))
        rolling = np.convolve(scores, np.ones(width), mode="valid")
        rows.append(
            {
                "class_name": class_name,
                "source": source,
                "sequence_id": identifier,
                "length": len(sequence),
                "prediction": prediction,
                "total_contribution": float(scores.sum()),
                "max_absolute_window_contribution": float(np.max(np.abs(rolling))),
            }
        )
    return pd.DataFrame(rows)


def parse_model_header(path: Path) -> dict[str, Any]:
    header: dict[str, Any] = {}
    with path.open() as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line == "SV":
                break
            key, *values = line.split()
            header[key] = values[0] if len(values) == 1 else values
    for key in ("L", "k", "d", "norc", "nr_class", "total_sv"):
        if key in header and not isinstance(header[key], list):
            header[key] = int(header[key])
    return header


def audit_support_vectors(path: Path) -> dict[str, Any]:
    in_support_vectors = False
    total = 0
    ambiguous = 0
    symbols: set[str] = set()
    with path.open() as handle:
        for line in handle:
            if not in_support_vectors:
                if line.rstrip("\n") == "SV":
                    in_support_vectors = True
                continue
            fields = line.split()
            if not fields:
                continue
            total += 1
            sequence = fields[-1].upper()
            unknown = set(sequence) - set(DNA)
            if unknown:
                ambiguous += 1
                symbols.update(unknown)
    return {
        "total": total,
        "ambiguous": ambiguous,
        "ambiguous_symbols": sorted(symbols),
        "gkmexplain_behavior": "non-ACGT support-vector bases are interpreted as A",
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_version(name: str) -> str:
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def summarize(args: argparse.Namespace) -> None:
    started = time.time()
    inputs = [args.model, args.positive_fasta, args.negative_fasta, args.motif_meme]
    for class_name in ("pos", "neg"):
        inputs.extend(
            [
                args.raw_dir / f"test_{class_name}_mode0.txt",
                args.raw_dir / f"test_{class_name}_mode1.txt",
                args.raw_dir / f"test_{class_name}_dinucleotide_shuffled.fa",
                args.raw_dir / f"test_{class_name}_null_mode0.txt",
            ]
        )
    for path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")
    if args.n_jobs <= 0 or args.window_size <= 0 or args.score_shards <= 0:
        raise ValueError("--n-jobs, --window-size, and --score-shards must be positive")
    if not 0 < args.target_seqlet_fdr <= 1 or not 0 < args.tomtom_pvalue <= 1:
        raise ValueError("FDR and p-value thresholds must be in (0, 1]")

    header = parse_model_header(args.model)
    kmer_length = args.kmer_length or int(header["L"])
    if kmer_length <= 0:
        raise ValueError("K-mer length must be positive")

    output_files = [
        "manifest.tsv", "kmer_importance.tsv", "motif_patterns.tsv",
        "null_diagnostics.tsv", "metadata.json", "positive_modisco.h5",
        "negative_modisco.h5",
    ]
    existing = [args.output_dir / name for name in output_files if (args.output_dir / name).exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Derived outputs already exist: {existing[0]}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".gkmexplain_summary.", dir=args.output_dir))
    try:
        motifs = read_meme(str(args.motif_meme))
        all_manifest: list[pd.DataFrame] = []
        all_pattern_tables: list[pd.DataFrame] = []
        all_diagnostics: list[pd.DataFrame] = []
        kmer_accumulators: dict[str, list[float]] = defaultdict(lambda: [0.0] * 10)

        class_specs = (
            ("positive", "pos", args.positive_fasta, 0),
            ("negative", "neg", args.negative_fasta, 1),
        )
        for class_name, short_name, fasta, class_index in class_specs:
            records = read_fasta(fasta)
            observed_path = args.raw_dir / f"test_{short_name}_mode0.txt"
            hypothetical_path = args.raw_dir / f"test_{short_name}_mode1.txt"
            X, H, manifest_rows, position_scores = load_class_arrays(
                records, observed_path, hypothetical_path, args.score_shards
            )
            manifest = pd.DataFrame(manifest_rows)
            manifest.insert(1, "class_name", class_name)
            all_manifest.append(manifest)
            accumulate_kmers(
                records, position_scores, class_index, kmer_length, kmer_accumulators
            )

            pos_patterns, neg_patterns = tfmodisco.TFMoDISco(
                one_hot=X,
                hypothetical_contribs=H,
                sliding_window_size=args.window_size,
                target_seqlet_fdr=args.target_seqlet_fdr,
                verbose=True,
            )
            modisco_io.save_hdf5(
                work / f"{class_name}_modisco.h5",
                pos_patterns,
                neg_patterns,
                window_size=args.window_size,
            )
            all_pattern_tables.append(
                annotate_patterns(
                    class_name,
                    pos_patterns,
                    neg_patterns,
                    motifs,
                    args.n_jobs,
                    args.tomtom_pvalue,
                )
            )

            all_diagnostics.append(
                window_diagnostics(
                    class_name,
                    "observed",
                    records,
                    observed_path,
                    args.window_size,
                    args.score_shards,
                )
            )
            null_fasta = args.raw_dir / f"test_{short_name}_dinucleotide_shuffled.fa"
            all_diagnostics.append(
                window_diagnostics(
                    class_name,
                    "dinucleotide_shuffled",
                    read_fasta(null_fasta),
                    args.raw_dir / f"test_{short_name}_null_mode0.txt",
                    args.window_size,
                    args.score_shards,
                )
            )
            del X, H, position_scores

        pd.concat(all_manifest, ignore_index=True).to_csv(
            work / "manifest.tsv", sep="\t", index=False
        )
        kmer_table(kmer_accumulators).to_csv(
            work / "kmer_importance.tsv", sep="\t", index=False
        )
        pd.concat(all_pattern_tables, ignore_index=True).to_csv(
            work / "motif_patterns.tsv", sep="\t", index=False
        )
        diagnostics = pd.concat(all_diagnostics, ignore_index=True)
        diagnostics.to_csv(work / "null_diagnostics.tsv", sep="\t", index=False)

        metadata = {
            "task": args.task,
            "model": str(args.model),
            "positive_fasta": str(args.positive_fasta),
            "negative_fasta": str(args.negative_fasta),
            "raw_dir": str(args.raw_dir),
            "motif_meme": str(args.motif_meme),
            "model_header": header,
            "support_vector_audit": audit_support_vectors(args.model),
            "model_compatibility": {
                "removed_header_from_temporary_copy": "norc 0" if header.get("norc") == 0 else None,
                "original_model_modified": False,
            },
            "parameters": {
                "kmer_length": kmer_length,
                "window_size": args.window_size,
                "target_seqlet_fdr": args.target_seqlet_fdr,
                "tomtom_pvalue": args.tomtom_pvalue,
                "null_sequences_per_class": 1000,
                "score_shards": args.score_shards,
            },
            "checksums": {
                "model": sha256(args.model),
                "positive_fasta": sha256(args.positive_fasta),
                "negative_fasta": sha256(args.negative_fasta),
                "motif_meme": sha256(args.motif_meme),
            },
            "software": {
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "modisco_lite": package_version("modisco"),
                "memelite": package_version("memelite"),
                "deeplift": package_version("deeplift"),
            },
            "slurm": {
                key: os.environ.get(key)
                for key in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID")
            },
            "elapsed_seconds": time.time() - started,
        }
        with (work / "metadata.json").open("w") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")

        for name in output_files:
            destination = args.output_dir / name
            if destination.exists():
                if destination.is_dir():
                    shutil.rmtree(destination)
                else:
                    destination.unlink()
            (work / name).replace(destination)
        print(json.dumps(metadata, indent=2, sort_keys=True))
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> None:
    args = parse_args()
    random.seed(getattr(args, "seed", 42))
    np.random.seed(getattr(args, "seed", 42))
    if args.command == "make-null":
        make_null(args)
    else:
        summarize(args)


if __name__ == "__main__":
    main()
