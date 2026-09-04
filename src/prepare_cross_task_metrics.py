#!/usr/bin/env python3
"""Combine native and cross-task EVN/SVN metrics for both model families."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from sklearn.metrics import average_precision_score, roc_auc_score

FRACTIONS = (10, 20, 40, 60, 80, 100)
TASKS = ("evn", "svn")


def count_fasta(path: Path) -> int:
    with path.open() as handle:
        return sum(line.startswith(">") for line in handle)


def metrics(labels: list[int], scores: list[float]) -> tuple[int, int, int, float, float]:
    return (len(labels), sum(labels), len(labels) - sum(labels),
            float(roc_auc_score(labels, scores)),
            float(average_precision_score(labels, scores)))


def read_native_lsgkm(path: Path, split: str) -> tuple[list[int], list[float]]:
    with path.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["split"] == split]
    return [int(row["label"]) for row in rows], [float(row["score"]) for row in rows]


def read_lsgkm_scores(pos: Path, neg: Path) -> tuple[list[int], list[float]]:
    def scores(path: Path) -> list[float]:
        with path.open() as handle:
            return [float(line.split()[-1]) for line in handle if line.strip()]
    positive = scores(pos)
    negative = scores(neg)
    return [1] * len(positive) + [0] * len(negative), positive + negative


def read_alphagenome(path: Path) -> tuple[list[int], list[float]]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [int(row["label"]) for row in rows], [float(row["logit"]) for row in rows]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path,
                        default=Path("results/cross_task/GRCh38/native_vs_cross_metrics.tsv"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = args.repository.resolve()
    rows: list[dict[str, object]] = []
    for model in ("lsgkm", "alphagenome"):
        model_dir = "ls-gkm" if model == "lsgkm" else "alphagenome"
        for model_task in TASKS:
            cross_task = "svn" if model_task == "evn" else "evn"
            for fraction in FRACTIONS:
                train_dir = root / "data/processed" / model_task / "GRCh38/learning_curves"
                n_train = count_fasta(train_dir / f"train_{fraction}_pos.fa") + count_fasta(train_dir / f"train_{fraction}_neg.fa")
                for evaluation_task, mode in ((model_task, "native"), (cross_task, "cross")):
                    for split in ("validation", "test"):
                        if model == "lsgkm" and mode == "native":
                            source = root / "results" / model_task / "GRCh38/ls-gkm/learning_curves" / f"{fraction}.tsv"
                            labels, scores = read_native_lsgkm(source, split)
                        elif model == "lsgkm":
                            source = root / "predictions/cross_task/GRCh38/ls-gkm" / f"{model_task}_model_on_{evaluation_task}" / str(fraction)
                            labels, scores = read_lsgkm_scores(source / f"{split}_pos_scores.txt", source / f"{split}_neg_scores.txt")
                        elif mode == "native":
                            source = root / "results" / model_task / "GRCh38/alphagenome/metrics.tsv"
                            with source.open() as handle:
                                candidates = [r for r in csv.DictReader(handle, delimiter="\t")
                                              if r["dataset"] == f"learning_curve_{fraction}" and r["split"] == split]
                            if len(candidates) != 1:
                                raise ValueError(f"Expected one native AlphaGenome row: {source}")
                            native = candidates[0]
                            rows.append(dict(model=model, model_task=model_task,
                                             evaluation_task=evaluation_task, evaluation_mode=mode,
                                             training_fraction=fraction, n_training_total=n_train,
                                             split=split, n_evaluation_total=int(native["n_evaluation_total"]),
                                             n_positive=int(native["n_evaluation_positive"]),
                                             n_negative=int(native["n_evaluation_negative"]),
                                             auroc=float(native["auroc"]), auprc=float(native["auprc"])))
                            continue
                        else:
                            source = root / "predictions/cross_task/GRCh38/alphagenome" / f"{model_task}_model_on_{evaluation_task}" / str(fraction) / f"{split}.tsv"
                            labels, scores = read_alphagenome(source)
                        n, pos, neg, auroc, auprc = metrics(labels, scores)
                        rows.append(dict(model=model, model_task=model_task,
                                         evaluation_task=evaluation_task, evaluation_mode=mode,
                                         training_fraction=fraction, n_training_total=n_train,
                                         split=split, n_evaluation_total=n,
                                         n_positive=pos, n_negative=neg,
                                         auroc=auroc, auprc=auprc))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
