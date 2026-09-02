#!/usr/bin/env python3
"""Validate and aggregate AlphaGenome EVN run artifacts for reporting."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable


FRACTIONS = (10, 20, 40, 60, 80, 100)
OUTERS = (1, 2, 3, 4)
METRIC_NAMES = ("loss", "accuracy", "auroc", "auprc", "precision", "recall", "f1")
CONFIG_FIELDS = (
    "sequence_length",
    "hidden_size",
    "dropout",
    "batch_size",
    "epochs",
    "learning_rate",
    "weight_decay",
    "patience",
    "seed",
    "threshold",
    "max_shift",
    "rc_prob",
    "shift_prob",
    "device_name",
    "torch_version",
    "numpy_version",
    "weights",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aggregate completed AlphaGenome learning-curve and outer runs."
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=Path("models/evn/GRCh38/alphagenome"),
        help="Directory containing learning_curves/ and outer1..outer4 (default: %(default)s).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/evn/GRCh38/alphagenome"),
        help="Destination for metrics.tsv and history.tsv (default: %(default)s).",
    )
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty artifact: {path}")
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty artifact: {path}")
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    return rows


def count_fasta(path_value: str) -> int:
    path = Path(path_value)
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty FASTA referenced by config: {path}")
    count = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    if count == 0:
        raise ValueError(f"No FASTA records found in {path}")
    return count


def require_metric_block(block: Any, path: Path, split: str) -> dict[str, float]:
    if not isinstance(block, dict):
        raise ValueError(f"Missing {split} metric object in {path}")
    missing = [name for name in METRIC_NAMES if name not in block]
    if missing:
        raise ValueError(f"Missing {split} metrics in {path}: {', '.join(missing)}")
    metrics = {name: float(block[name]) for name in METRIC_NAMES}
    if not all(math.isfinite(value) for value in metrics.values()):
        raise ValueError(f"Non-finite {split} metric in {path}")
    return metrics


def run_specs(model_dir: Path) -> Iterable[tuple[str, str, int | None, Path]]:
    for fraction in FRACTIONS:
        dataset = f"learning_curve_{fraction}"
        yield dataset, "learning_curve", fraction, model_dir / "learning_curves" / str(fraction) / dataset
    for outer in OUTERS:
        dataset = f"outer{outer}"
        yield dataset, "outer", None, model_dir / dataset / dataset


def atomic_write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main() -> None:
    args = parse_args()
    metrics_rows: list[dict[str, Any]] = []
    history_rows: list[dict[str, Any]] = []
    reference_config: dict[str, Any] | None = None

    for dataset, dataset_type, fraction, prefix in run_specs(args.model_dir):
        config_path = Path(str(prefix) + ".config.json")
        history_path = Path(str(prefix) + ".history.tsv")
        validation_path = Path(str(prefix) + ".validation_metrics.json")
        test_path = Path(str(prefix) + ".test_metrics.json")
        checkpoint_path = Path(str(prefix) + ".best_head.pt")
        if not checkpoint_path.is_file() or checkpoint_path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty checkpoint: {checkpoint_path}")

        config = read_json(config_path)
        validation = read_json(validation_path)
        test = read_json(test_path)
        history = read_tsv(history_path)

        missing_config = [field for field in CONFIG_FIELDS if field not in config]
        if missing_config:
            raise ValueError(f"Missing config values in {config_path}: {', '.join(missing_config)}")
        selected_config = {field: config[field] for field in CONFIG_FIELDS}
        if reference_config is None:
            reference_config = selected_config
        elif selected_config != reference_config:
            differing = [field for field in CONFIG_FIELDS if selected_config[field] != reference_config[field]]
            raise ValueError(f"Frozen configuration mismatch in {config_path}: {', '.join(differing)}")

        best_epoch = int(validation.get("best_epoch", -1))
        if validation.get("selection_metric") != "validation_auprc":
            raise ValueError(f"Unexpected checkpoint selection metric in {validation_path}")
        if int(test.get("best_epoch", -2)) != best_epoch:
            raise ValueError(f"Best-epoch mismatch between {validation_path} and {test_path}")
        if test.get("checkpoint_selection") != "validation_auprc":
            raise ValueError(f"Unexpected test checkpoint selection in {test_path}")

        validation_metrics = require_metric_block(validation.get("validation"), validation_path, "validation")
        test_metrics = require_metric_block(test.get("test"), test_path, "test")
        elapsed_seconds = float(test.get("elapsed_seconds", float("nan")))
        if not math.isfinite(elapsed_seconds) or elapsed_seconds <= 0:
            raise ValueError(f"Invalid elapsed_seconds in {test_path}")

        history_by_epoch = {int(row["epoch"]): row for row in history}
        if best_epoch not in history_by_epoch:
            raise ValueError(f"Best epoch {best_epoch} absent from {history_path}")
        recorded_best = float(history_by_epoch[best_epoch]["validation_auprc"])
        if not math.isclose(recorded_best, validation_metrics["auprc"], rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Best validation AUPRC does not match history in {history_path}")

        train_positive = count_fasta(str(config["train_pos"]))
        train_negative = count_fasta(str(config["train_neg"]))
        validation_positive = count_fasta(str(config["validation_pos"]))
        validation_negative = count_fasta(str(config["validation_neg"]))
        test_positive = count_fasta(str(config["test_pos"]))
        test_negative = count_fasta(str(config["test_neg"]))
        if train_positive != train_negative:
            raise ValueError(f"Training class imbalance in {dataset}")
        if validation_positive != validation_negative or test_positive != test_negative:
            raise ValueError(f"Evaluation class imbalance in {dataset}")

        common: dict[str, Any] = {
            "model": "alphagenome",
            "task": "GRCh38",
            "dataset": dataset,
            "dataset_type": dataset_type,
            "training_fraction": "" if fraction is None else fraction,
            "n_training_positive": train_positive,
            "n_training_negative": train_negative,
            "n_training_total": train_positive + train_negative,
            "best_epoch": best_epoch,
            "epochs_completed": len(history),
            "selection_metric": "validation_auprc",
            "elapsed_seconds": elapsed_seconds,
            **selected_config,
        }
        for split, values, positive, negative in (
            ("validation", validation_metrics, validation_positive, validation_negative),
            ("test", test_metrics, test_positive, test_negative),
        ):
            metrics_rows.append(
                {
                    **common,
                    "split": split,
                    "n_evaluation_positive": positive,
                    "n_evaluation_negative": negative,
                    "n_evaluation_total": positive + negative,
                    **values,
                }
            )

        for row in history:
            required = [
                "epoch",
                *(f"{split}_{metric}" for split in ("train", "validation") for metric in METRIC_NAMES),
            ]
            missing = [field for field in required if field not in row]
            if missing:
                raise ValueError(f"Missing history columns in {history_path}: {', '.join(missing)}")
            history_rows.append(
                {
                    "model": "alphagenome",
                    "task": "GRCh38",
                    "dataset": dataset,
                    "dataset_type": dataset_type,
                    "training_fraction": "" if fraction is None else fraction,
                    "n_training_total": train_positive + train_negative,
                    "best_epoch": best_epoch,
                    **{field: row[field] for field in required},
                }
            )

    metric_fields = list(metrics_rows[0])
    history_fields = list(history_rows[0])
    metrics_path = args.output_dir / "metrics.tsv"
    history_output_path = args.output_dir / "history.tsv"
    atomic_write_tsv(metrics_path, metrics_rows, metric_fields)
    atomic_write_tsv(history_output_path, history_rows, history_fields)
    print(f"Wrote {len(metrics_rows)} rows to {metrics_path}")
    print(f"Wrote {len(history_rows)} rows to {history_output_path}")
    print("AlphaGenome artifact validation: PASS")


if __name__ == "__main__":
    main()
