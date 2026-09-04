#!/usr/bin/env python3
"""Compare EVN/SVN AlphaGenome head connection strength in encoder space.

For input feature j, importance is sum_h |w2_h * w1_hj| multiplied by the
absolute LayerNorm scale for the corresponding encoder channel. This is a
model-global path-strength measure. It deliberately does not claim to capture
sample-dependent ReLU gating or nucleotide-level causal attribution.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import torch

FRACTIONS = (10, 20, 40, 60, 80, 100)
TASKS = ("evn", "svn")
ENCODER_DIM = 1536


def rankdata(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    return ranks


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


def importance(checkpoint_path: Path) -> np.ndarray:
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    state = checkpoint["head"]
    w1 = state["fc1.weight"].float().numpy()
    w2 = state["fc2.weight"].float().numpy().reshape(-1)
    gamma = state["norm.weight"].float().numpy()
    positions = int(checkpoint["encoder_positions"])
    if w1.shape[1] != positions * ENCODER_DIM or gamma.shape != (ENCODER_DIM,):
        raise ValueError(f"Unexpected head dimensions in {checkpoint_path}")
    path_strength = (np.abs(w2)[:, None] * np.abs(w1)).sum(axis=0)
    path_strength *= np.tile(np.abs(gamma), positions)
    return path_strength


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-root", type=Path, default=Path("models"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("results/cross_task/GRCh38/alphagenome_features"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    feature_rows = []
    similarity_rows = []
    for fraction in FRACTIONS:
        vectors = {}
        for task in TASKS:
            checkpoint = (args.models_root / task / "GRCh38/alphagenome/learning_curves" /
                          str(fraction) / f"learning_curve_{fraction}.best_head.pt")
            vector = importance(checkpoint)
            vectors[task] = vector
            normalized = vector / vector.sum()
            for index, value in enumerate(vector):
                feature_rows.append({
                    "task": task, "training_fraction": fraction,
                    "feature_index": index,
                    "encoder_position": index // ENCODER_DIM,
                    "encoder_channel": index % ENCODER_DIM,
                    "path_strength": float(value),
                    "normalized_importance": float(normalized[index]),
                })
        a, b = vectors["evn"], vectors["svn"]
        row = {"training_fraction": fraction,
               "cosine_similarity": cosine(a, b),
               "spearman_correlation": spearman(a, b)}
        for percent in (1, 5, 10):
            n_top = max(1, len(a) * percent // 100)
            top_a = set(np.argpartition(a, -n_top)[-n_top:])
            top_b = set(np.argpartition(b, -n_top)[-n_top:])
            overlap = len(top_a & top_b)
            row[f"top_{percent}_overlap"] = overlap
            row[f"top_{percent}_overlap_fraction"] = overlap / n_top
        similarity_rows.append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for filename, rows in (("feature_importance.tsv", feature_rows),
                           ("feature_similarity.tsv", similarity_rows)):
        path = args.output_dir / filename
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {path}")


if __name__ == "__main__":
    main()
