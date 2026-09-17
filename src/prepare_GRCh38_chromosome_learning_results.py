#!/usr/bin/env python3
"""Aggregate cumulative chromosome-learning metrics for LS-GKM or AlphaGenome."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

CHROMS = tuple(str(i) for i in range(3, 23)) + ("X",)
METRICS = ("loss", "accuracy", "auroc", "auprc", "precision", "recall", "f1")


def fasta_count(path: Path) -> int:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty FASTA: {path}")
    with path.open(encoding="utf-8") as handle:
        count = sum(line.startswith(">") for line in handle)
    if not count:
        raise ValueError(f"No FASTA records in {path}")
    return count


def read_json(path: Path) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty JSON: {path}")
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def read_scores(path: Path) -> list[float]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty score file: {path}")
    values = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed score row")
            try:
                values.append(float(fields[-1]))
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid score") from error
    return values


def binary_metrics(positive: list[float], negative: list[float]) -> dict[str, float]:
    scores = positive + negative
    labels = [1] * len(positive) + [0] * len(negative)
    if not positive or not negative:
        raise ValueError("Both classes must be non-empty")
    order = sorted(range(len(scores)), key=lambda i: scores[i])
    ranks = [0.0] * len(scores)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and scores[order[end]] == scores[order[start]]:
            end += 1
        average_rank = (start + 1 + end) / 2
        for position in range(start, end):
            ranks[order[position]] = average_rank
        start = end
    n_pos, n_neg = len(positive), len(negative)
    rank_sum = sum(rank for rank, label in zip(ranks, labels) if label == 1)
    auroc = (rank_sum - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
    ranked = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    true_positive = 0
    precision_sum = 0.0
    for rank, index in enumerate(ranked, 1):
        if labels[index] == 1:
            true_positive += 1
            precision_sum += true_positive / rank
    auprc = precision_sum / n_pos
    predictions = [int(score >= 0) for score in scores]
    tp = sum(p == 1 and y == 1 for p, y in zip(predictions, labels))
    tn = sum(p == 0 and y == 0 for p, y in zip(predictions, labels))
    fp = sum(p == 1 and y == 0 for p, y in zip(predictions, labels))
    fn = sum(p == 0 and y == 1 for p, y in zip(predictions, labels))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    return {
        "loss": math.nan,
        "accuracy": (tp + tn) / len(labels),
        "auroc": auroc,
        "auprc": auprc,
        "precision": precision,
        "recall": recall,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", choices=("evn", "svn"), required=True)
    parser.add_argument("--model", choices=("lsgkm", "alphagenome"), required=True)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    repo = args.repo.resolve()
    training = repo / f"data/processed/{args.dataset}/GRCh38/training"
    model_name = "ls-gkm" if args.model == "lsgkm" else "alphagenome"
    output = args.output or (
        repo / f"results/{args.dataset}/GRCh38/{model_name}/"
        "chromosome_learning_curves/metrics.tsv"
    )
    evaluation_counts = {
        split: {cls: fasta_count(training / f"{split}_{cls}.fa") for cls in ("pos", "neg")}
        for split in ("train_100", "validation", "test")
    }
    cumulative_pos = cumulative_neg = 0
    rows: list[dict[str, object]] = []
    for checkpoint_index, chrom in enumerate(CHROMS, 1):
        cumulative_pos += fasta_count(training / f"chromosomes/chrom{chrom}/train_pos.fa")
        cumulative_neg += fasta_count(training / f"chromosomes/chrom{chrom}/train_neg.fa")
        if cumulative_pos != cumulative_neg:
            raise ValueError(f"Unbalanced cumulative training data at chr{chrom}")
        name = f"through_chr{chrom}"
        if args.model == "alphagenome":
            prefix = repo / f"models/{args.dataset}/GRCh38/alphagenome/chromosome_learning_curves/{name}/{name}"
            config = read_json(Path(str(prefix) + ".config.json"))
            validation_json = read_json(Path(str(prefix) + ".validation_metrics.json"))
            best_epoch = validation_json["best_epoch"]
            artifacts = {
                "train_100": read_json(Path(str(prefix) + ".train_100_metrics.json"))["metrics"],
                "validation": validation_json["validation"],
                "test": read_json(Path(str(prefix) + ".test_metrics.json"))["test"],
            }
        else:
            prediction_dir = repo / f"predictions/{args.dataset}/GRCh38/ls-gkm/chromosome_learning_curves/{name}"
            config = {}
            best_epoch = ""
            artifacts = {
                split: binary_metrics(
                    read_scores(prediction_dir / f"{split}_pos_scores.txt"),
                    read_scores(prediction_dir / f"{split}_neg_scores.txt"),
                )
                for split in ("train_100", "validation", "test")
            }
        for split, metrics in artifacts.items():
            n_pos = evaluation_counts[split]["pos"]
            n_neg = evaluation_counts[split]["neg"]
            row: dict[str, object] = {
                "model": args.model, "task": "GRCh38", "dataset": args.dataset,
                "negative_kind": "genomic", "checkpoint": name,
                "checkpoint_index": checkpoint_index, "endpoint_chromosome": f"chr{chrom}",
                "training_chromosomes": ",".join(f"chr{x}" for x in CHROMS[:checkpoint_index]),
                "n_training_positive": cumulative_pos, "n_training_negative": cumulative_neg,
                "n_training_total": cumulative_pos + cumulative_neg, "split": split,
                "n_evaluation_positive": n_pos, "n_evaluation_negative": n_neg,
                "n_evaluation_total": n_pos + n_neg, "best_epoch": best_epoch,
            }
            row.update({metric: metrics.get(metric, math.nan) for metric in METRICS})
            for field in ("sequence_length", "hidden_size", "dropout", "batch_size",
                          "epochs", "learning_rate", "patience", "seed", "weights"):
                row[field] = config.get(field, "")
            rows.append(row)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(output) + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(output)
    print(f"Wrote {len(rows)} rows to {output}")
    print("Chromosome-learning metric aggregation: PASS")


if __name__ == "__main__":
    main()
