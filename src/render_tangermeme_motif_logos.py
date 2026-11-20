#!/usr/bin/env python3
"""Render paired model-derived and HOCOMOCO logos for the Tangermeme report."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from memelite.io import read_meme
from tangermeme.plot import plot_logo, plot_pwm

from tangermeme_analysis import sample_and_encode


TASKS = ("evn", "svn", "evs")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-branch", required=True)
    parser.add_argument("--motif-meme", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-per-task", type=int, default=5)
    return parser.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def reverse_complement_matrix(matrix: np.ndarray) -> np.ndarray:
    return matrix[[3, 2, 1, 0], ::-1]


def align_contributions(
    X: torch.Tensor,
    attribution: torch.Tensor,
    sites: pd.DataFrame,
    pwm: np.ndarray,
) -> np.ndarray:
    """Align candidate contributions to the best PWM strand and offset."""
    pwm = np.asarray(pwm, dtype=np.float64)
    motif_length = pwm.shape[1]
    log_pwm = np.log2(np.clip(pwm, 1e-7, None) / 0.25)
    reverse_log_pwm = reverse_complement_matrix(log_pwm)
    aligned: list[np.ndarray] = []

    for site in sites.itertuples(index=False):
        example_idx = int(site.example_idx)
        sequence = X[example_idx].numpy().astype(np.float64, copy=False)
        contribution = attribution[example_idx].numpy().astype(np.float64, copy=False)
        sequence_length = sequence.shape[1]
        first = max(0, int(site.start) - motif_length + 1)
        last = min(sequence_length - motif_length, int(site.end) - 1)
        if last < first:
            continue

        best_score = -np.inf
        best_start = first
        best_reverse = False
        for start in range(first, last + 1):
            window = sequence[:, start : start + motif_length]
            forward_score = float((window * log_pwm).sum())
            reverse_score = float((window * reverse_log_pwm).sum())
            if forward_score > best_score:
                best_score = forward_score
                best_start = start
                best_reverse = False
            if reverse_score > best_score:
                best_score = reverse_score
                best_start = start
                best_reverse = True

        window = contribution[:, best_start : best_start + motif_length]
        if best_reverse:
            window = reverse_complement_matrix(window)
        aligned.append(window)

    if not aligned:
        raise RuntimeError("No candidate sites could be aligned to the motif")
    return np.stack(aligned).mean(axis=0)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    motifs = read_meme(args.motif_meme)
    records: list[dict[str, object]] = []

    for task in TASKS:
        result_dir = (
            Path("results")
            / task
            / "GRCh38/alphagenome/real"
            / args.result_branch
            / "through_chrX"
        )
        summary = pd.read_csv(result_dir / "motif_cluster_summary.tsv", sep="\t")
        candidates = pd.read_csv(result_dir / "candidate_region_tests.tsv", sep="\t")
        metadata = pd.read_json(result_dir / "metadata.json", typ="series")
        dataset = metadata["dataset"]
        parameters = metadata["parameters"]
        X, manifest = sample_and_encode(
            Path(metadata["positive_fasta"]),
            Path(metadata["negative_fasta"]),
            int(parameters["n_per_class"]),
            int(dataset["sequence_length"]),
            int(parameters["seed"]),
            all_sequences=bool(dataset["all_sequences"]),
        )
        stored_manifest = pd.read_csv(result_dir / "manifest.tsv", sep="\t")
        if not manifest["sequence_id"].equals(stored_manifest["sequence_id"]):
            raise RuntimeError(f"{task}: reconstructed sequence order differs from manifest")
        attribution = torch.load(
            result_dir / "attributions.pt", map_location="cpu", weights_only=True
        )
        eligible = summary.loc[summary["n_sites"] >= 2].copy()
        eligible = eligible.sort_values(
            ["n_important_regions", "mean_ablation_delta"],
            ascending=[False, False],
            key=lambda column: column.abs()
            if column.name == "mean_ablation_delta"
            else column,
        )
        selected = eligible.head(args.top_per_task)
        enriched = eligible.loc[eligible["enriched_cluster"].astype(str).str.lower().isin(("true", "1"))]
        selected = pd.concat([selected, enriched]).drop_duplicates(
            ["motif", "sign"], keep="first"
        )

        for rank, row in enumerate(selected.itertuples(index=False), start=1):
            motif_name = str(row.motif)
            if motif_name not in motifs:
                raise KeyError(f"{motif_name!r} is absent from {args.motif_meme}")
            pwm = motifs[motif_name]
            sites = candidates.loc[
                (candidates["motif"] == motif_name)
                & (candidates["sign"] == row.sign)
            ]
            mean_contribution = align_contributions(X, attribution, sites, pwm)
            width = max(8.8, 0.82 * pwm.shape[1])
            fig, axes = plt.subplots(1, 2, figsize=(width, 2.7))
            plot_logo(mean_contribution, ax=axes[0])
            axes[0].set_title(
                f"Model-derived mean contribution (n={len(sites)})"
            )
            axes[0].set_xlabel("Aligned motif position")
            axes[0].set_ylabel("Mean attribution")
            plot_pwm(pwm, ax=axes[1], name=str(row.representative_gene))
            axes[1].set_title(f"HOCOMOCO reference: {row.representative_gene}")
            axes[1].set_xlabel("Motif position")
            axes[1].set_ylabel("Information (bits)")
            fig.tight_layout()
            output = args.output_dir / (
                f"{task}_{rank:02d}_{safe_name(str(row.representative_gene))}_"
                f"{safe_name(str(row.sign))}.png"
            )
            fig.savefig(output, dpi=180, bbox_inches="tight")
            plt.close(fig)
            records.append(
                {
                    "task": task.upper(),
                    "rank": rank,
                    "motif": motif_name,
                    "representative_gene": row.representative_gene,
                    "sign": row.sign,
                    "n_sites": int(row.n_sites),
                    "n_important_regions": int(row.n_important_regions),
                    "enriched_cluster": bool(row.enriched_cluster),
                    "n_aligned_sites": len(sites),
                    "png": str(output.resolve()),
                }
            )

    pd.DataFrame(records).to_csv(
        args.output_dir / "motif_logo_manifest.tsv", sep="\t", index=False
    )


if __name__ == "__main__":
    main()
