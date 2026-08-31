#!/usr/bin/env python3
"""Train one frozen-AlphaGenome enhancer-vs-negative classifier.

The command deliberately resembles ``gkmtrain``: it receives explicit positive
and negative training FASTAs plus an output prefix. Two validation FASTAs are
also required because validation AUPRC selects the checkpoint and controls early
stopping. Two test FASTAs are accepted for the final assessment, but they are not
loaded until training has ended and the validation-selected checkpoint is restored.

Example
-------
src/train_alphagenome.py \
    data/processed/evn/GRCh38/learning_curves/train_10_pos.fa \
    data/processed/evn/GRCh38/learning_curves/train_10_neg.fa \
    data/processed/evn/GRCh38/training/validation_pos.fa \
    data/processed/evn/GRCh38/training/validation_neg.fa \
    data/processed/evn/GRCh38/training/test_pos.fa \
    data/processed/evn/GRCh38/training/test_neg.fa \
    models/evn/GRCh38/alphagenome/learning_curves/10/learning_curve_10

The output prefix above produces:

* ``learning_curve_10.best_head.pt``
* ``learning_curve_10.config.json``
* ``learning_curve_10.history.tsv``
* ``learning_curve_10.validation_metrics.json``
* ``learning_curve_10.test_metrics.json``
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import random
import sys
import time
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from torch import Tensor
from torch.utils.data import DataLoader, Dataset
from tqdm import tqdm

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.training import create_lr_scheduler
from alphagenome_pytorch.extensions.finetuning.transfer import load_trunk, remove_all_heads
from alphagenome_pytorch.extensions.finetuning.utils import sequence_to_onehot


ENCODER_DIM = 1536
TRACKED_METRICS = ("loss", "accuracy", "auroc", "auprc", "precision", "recall", "f1")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train one binary enhancer classifier on a frozen pretrained "
            "AlphaGenome encoder."
        )
    )
    parser.add_argument("train_pos", type=Path, help="Positive training FASTA.")
    parser.add_argument("train_neg", type=Path, help="Negative training FASTA.")
    parser.add_argument("validation_pos", type=Path, help="Positive validation FASTA.")
    parser.add_argument("validation_neg", type=Path, help="Negative validation FASTA.")
    parser.add_argument("test_pos", type=Path, help="Positive held-out test FASTA.")
    parser.add_argument("test_neg", type=Path, help="Negative held-out test FASTA.")
    parser.add_argument(
        "output_prefix",
        type=Path,
        help="Output prefix, analogous to the final gkmtrain positional argument.",
    )
    parser.add_argument(
        "--weights",
        type=Path,
        default=Path("alphagenome_trial/finetune_files/weights/fold_0_weights.safetensors"),
        help="Pretrained AlphaGenome weights (default: %(default)s).",
    )
    parser.add_argument("--sequence-length", type=int, default=256)
    parser.add_argument("--hidden-size", type=int, default=1024)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--num-workers", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=10)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--weight-decay", type=float, default=0.0)
    parser.add_argument("--patience", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--max-shift", type=int, default=15)
    parser.add_argument("--rc-prob", type=float, default=0.5)
    parser.add_argument("--shift-prob", type=float, default=0.5)
    parser.add_argument("--no-amp", action="store_true", help="Disable bfloat16 autocast.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace artifacts that already exist for this output prefix.",
    )
    return parser.parse_args()


def output_paths(prefix: Path) -> dict[str, Path]:
    stem = str(prefix)
    return {
        "checkpoint": Path(stem + ".best_head.pt"),
        "config": Path(stem + ".config.json"),
        "history": Path(stem + ".history.tsv"),
        "validation_metrics": Path(stem + ".validation_metrics.json"),
        "test_metrics": Path(stem + ".test_metrics.json"),
    }


def validate_args(args: argparse.Namespace, outputs: dict[str, Path]) -> None:
    required = (
        args.train_pos,
        args.train_neg,
        args.validation_pos,
        args.validation_neg,
        args.test_pos,
        args.test_neg,
        args.weights,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required file(s):\n  " + "\n  ".join(missing))
    empty = [str(path) for path in required if path.stat().st_size == 0]
    if empty:
        raise ValueError("Required file(s) are empty:\n  " + "\n  ".join(empty))
    if args.sequence_length <= 0 or args.sequence_length % 128:
        raise ValueError("--sequence-length must be positive and divisible by 128")
    if args.hidden_size <= 0 or args.batch_size <= 0 or args.num_workers < 0:
        raise ValueError("hidden size and batch size must be positive; workers cannot be negative")
    if args.epochs <= 0 or args.patience <= 0:
        raise ValueError("--epochs and --patience must be positive")
    if args.learning_rate <= 0 or args.weight_decay < 0:
        raise ValueError("learning rate must be positive and weight decay cannot be negative")
    if not 0 <= args.dropout < 1:
        raise ValueError("--dropout must be in [0, 1)")
    if not 0 <= args.rc_prob <= 1 or not 0 <= args.shift_prob <= 1:
        raise ValueError("augmentation probabilities must be in [0, 1]")
    if args.max_shift < 0:
        raise ValueError("--max-shift cannot be negative")
    if not 0 <= args.threshold <= 1:
        raise ValueError("--threshold must be in [0, 1]")
    existing = [str(path) for path in outputs.values() if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Output artifact(s) already exist; use --overwrite to replace them:\n  "
            + "\n  ".join(existing)
        )


def require_cuda() -> None:
    """Print the PyTorch/GPU environment and fail before loading any data."""
    print(f"PyTorch: {torch.__version__}", flush=True)
    print(f"CUDA available: {torch.cuda.is_available()}", flush=True)
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable; submit this trainer to a GPU compute node")
    print(f"GPU: {torch.cuda.get_device_name(0)}", flush=True)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def seed_worker(worker_id: int) -> None:
    del worker_id
    worker_seed = torch.initial_seed() % (2**32)
    random.seed(worker_seed)
    np.random.seed(worker_seed)


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    sequence_parts: list[str] = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(sequence_parts).upper()))
                identifier = line[1:].split()[0]
                if not identifier:
                    raise ValueError(f"{path}:{line_number}: empty FASTA identifier")
                sequence_parts = []
            else:
                if identifier is None:
                    raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
                sequence_parts.append(line)
    if identifier is not None:
        records.append((identifier, "".join(sequence_parts).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    empty_ids = [identifier for identifier, sequence in records if not sequence]
    if empty_ids:
        raise ValueError(f"Empty sequence in {path}: {empty_ids[0]}")
    ids = [identifier for identifier, _ in records]
    if len(ids) != len(set(ids)):
        raise ValueError(f"Duplicate FASTA identifiers found in {path}")
    return records


class EnhancerDataset(Dataset):
    """Return center-adjusted one-hot sequences with enhancer=1 and negative=0."""

    def __init__(
        self,
        positive_fasta: Path,
        negative_fasta: Path,
        sequence_length: int,
        augment: bool = False,
        rc_prob: float = 0.5,
        shift_prob: float = 0.5,
        max_shift: int = 15,
    ) -> None:
        positives = read_fasta(positive_fasta)
        negatives = read_fasta(negative_fasta)
        self.sequences = [sequence for _, sequence in positives + negatives]
        self.labels = np.concatenate(
            (
                np.ones(len(positives), dtype=np.float32),
                np.zeros(len(negatives), dtype=np.float32),
            )
        )
        self.n_positive = len(positives)
        self.n_negative = len(negatives)
        if self.n_positive != self.n_negative:
            raise ValueError(
                f"Class imbalance: {self.n_positive:,} positives and "
                f"{self.n_negative:,} negatives"
            )
        self.sequence_length = sequence_length
        self.augment = augment
        self.rc_prob = rc_prob
        self.shift_prob = shift_prob
        self.max_shift = max_shift

    def __len__(self) -> int:
        return len(self.sequences)

    def __getitem__(self, index: int) -> tuple[Tensor, Tensor]:
        onehot = sequence_to_onehot(self.sequences[index]).astype(np.float32)
        onehot = self._center_pad_or_crop(onehot)
        if self.augment and torch.rand(()).item() < self.shift_prob:
            shift = int(torch.randint(-self.max_shift, self.max_shift + 1, ()).item())
            onehot = self._shift_with_zeros(onehot, shift)
        if self.augment and torch.rand(()).item() < self.rc_prob:
            onehot = onehot[::-1, :][:, [3, 2, 1, 0]].copy()
        return torch.from_numpy(onehot), torch.tensor(self.labels[index])

    def _center_pad_or_crop(self, onehot: np.ndarray) -> np.ndarray:
        length = onehot.shape[0]
        if length < self.sequence_length:
            difference = self.sequence_length - length
            left = difference // 2
            return np.pad(onehot, ((left, difference - left), (0, 0)))
        if length > self.sequence_length:
            start = (length - self.sequence_length) // 2
            return onehot[start : start + self.sequence_length]
        return onehot

    @staticmethod
    def _shift_with_zeros(onehot: np.ndarray, shift: int) -> np.ndarray:
        shifted = np.zeros_like(onehot)
        if shift > 0:
            shifted[shift:] = onehot[:-shift]
        elif shift < 0:
            shifted[:shift] = onehot[-shift:]
        else:
            shifted[:] = onehot
        return shifted


class EnhancerHead(nn.Module):
    """MPRA-style MLP returning one unbounded binary-classification logit."""

    def __init__(self, n_positions: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        self.norm = nn.LayerNorm(ENCODER_DIM)
        self.fc1 = nn.Linear(n_positions * ENCODER_DIM, hidden_size)
        self.dropout = nn.Dropout(dropout)
        self.fc2 = nn.Linear(hidden_size, 1)

    def forward(self, encoder_output: Tensor) -> Tensor:
        x = self.norm(encoder_output)
        x = x.flatten(1)
        x = F.relu(self.fc1(x))
        x = self.dropout(x)
        return self.fc2(x).squeeze(-1)


def amp_context(device: torch.device, enabled: bool):
    if enabled and device.type == "cuda":
        return torch.autocast(device_type="cuda", dtype=torch.bfloat16)
    return nullcontext()


def normalize_encoder_output(output: Tensor) -> Tensor:
    """Normalize encoder output to (batch, positions, 1536)."""
    if output.ndim != 3:
        raise ValueError(f"Expected 3D encoder output, received {tuple(output.shape)}")
    if output.shape[-1] == ENCODER_DIM:
        return output
    if output.shape[1] == ENCODER_DIM:
        return output.transpose(1, 2)
    raise ValueError(f"Cannot identify encoder dimension in shape {tuple(output.shape)}")


def encode(model: nn.Module, sequences: Tensor, device: torch.device, use_amp: bool) -> Tensor:
    organism = torch.zeros(sequences.shape[0], dtype=torch.long, device=device)
    with torch.no_grad(), amp_context(device, use_amp):
        output = model(sequences, organism, encoder_only=True)["encoder_output"]
    return normalize_encoder_output(output)


def calculate_metrics(
    probabilities: Tensor,
    targets: Tensor,
    average_loss: float,
    threshold: float,
) -> dict[str, float]:
    scores = probabilities.numpy()
    truth = targets.numpy().astype(int)
    predicted = (scores >= threshold).astype(int)
    if np.unique(truth).size != 2:
        raise ValueError("AUROC and AUPRC require both classes")
    return {
        "loss": float(average_loss),
        "accuracy": float(accuracy_score(truth, predicted)),
        "auroc": float(roc_auc_score(truth, scores)),
        "auprc": float(average_precision_score(truth, scores)),
        "precision": float(precision_score(truth, predicted, zero_division=0)),
        "recall": float(recall_score(truth, predicted, zero_division=0)),
        "f1": float(f1_score(truth, predicted, zero_division=0)),
    }


def train_epoch(
    model: nn.Module,
    head: EnhancerHead,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    scheduler,
    device: torch.device,
    use_amp: bool,
    threshold: float,
) -> dict[str, float]:
    model.eval()
    head.train()
    total_loss = 0.0
    total_examples = 0
    probabilities: list[Tensor] = []
    targets_all: list[Tensor] = []
    for sequences, targets in tqdm(loader, desc="train", leave=False):
        sequences = sequences.to(device, non_blocking=True)
        targets = targets.to(device, non_blocking=True).float()
        encoder_output = encode(model, sequences, device, use_amp)
        optimizer.zero_grad(set_to_none=True)
        with amp_context(device, use_amp):
            logits = head(encoder_output.detach())
            loss = F.binary_cross_entropy_with_logits(logits.float(), targets)
        loss.backward()
        optimizer.step()
        scheduler.step()
        batch_size = targets.numel()
        total_loss += loss.item() * batch_size
        total_examples += batch_size
        probabilities.append(torch.sigmoid(logits.detach().float()).cpu())
        targets_all.append(targets.cpu())
    return calculate_metrics(
        torch.cat(probabilities),
        torch.cat(targets_all),
        total_loss / total_examples,
        threshold,
    )


@torch.no_grad()
def evaluate(
    model: nn.Module,
    head: EnhancerHead,
    loader: DataLoader,
    device: torch.device,
    use_amp: bool,
    threshold: float,
    description: str,
) -> dict[str, float]:
    model.eval()
    head.eval()
    total_loss = 0.0
    total_examples = 0
    probabilities: list[Tensor] = []
    targets_all: list[Tensor] = []
    for sequences, targets in tqdm(loader, desc=description, leave=False):
        sequences = sequences.to(device, non_blocking=True)
        targets = targets.to(device, non_blocking=True).float()
        encoder_output = encode(model, sequences, device, use_amp)
        with amp_context(device, use_amp):
            logits = head(encoder_output)
            loss = F.binary_cross_entropy_with_logits(logits.float(), targets)
        batch_size = targets.numel()
        total_loss += loss.item() * batch_size
        total_examples += batch_size
        probabilities.append(torch.sigmoid(logits.float()).cpu())
        targets_all.append(targets.cpu())
    return calculate_metrics(
        torch.cat(probabilities),
        torch.cat(targets_all),
        total_loss / total_examples,
        threshold,
    )


def build_loader(
    dataset: EnhancerDataset,
    batch_size: int,
    workers: int,
    shuffle: bool,
    device: torch.device,
    generator: torch.Generator,
) -> DataLoader:
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=workers,
        pin_memory=device.type == "cuda",
        worker_init_fn=seed_worker,
        generator=generator if shuffle else None,
    )


def write_json(path: Path, value: object) -> None:
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def write_history(path: Path, history: list[dict[str, float]]) -> None:
    fields = ["epoch"] + [
        f"{split}_{metric}"
        for split in ("train", "validation")
        for metric in TRACKED_METRICS
    ]
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(history)
    temporary.replace(path)


def save_checkpoint(path: Path, payload: dict[str, object]) -> None:
    temporary = Path(str(path) + ".tmp")
    torch.save(payload, temporary)
    temporary.replace(path)


def configuration(args: argparse.Namespace) -> dict[str, object]:
    config = {
        key: str(value) if isinstance(value, Path) else value
        for key, value in vars(args).items()
    }
    config.update(
        {
            "torch_version": str(torch.__version__),
            "numpy_version": np.__version__,
            "python_version": sys.version,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "slurm_array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
            "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        }
    )
    return config


def main() -> None:
    args = parse_args()
    require_cuda()
    outputs = output_paths(args.output_prefix)
    validate_args(args, outputs)
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    seed_everything(args.seed)
    device = torch.device("cuda")
    use_amp = not args.no_amp
    started = time.time()

    config = configuration(args)
    config["device_name"] = torch.cuda.get_device_name(0) if device.type == "cuda" else "cpu"
    config["output_artifacts"] = {key: str(path) for key, path in outputs.items()}
    write_json(outputs["config"], config)
    print(json.dumps(config, indent=2, sort_keys=True), flush=True)

    print("Loading FASTA datasets...", flush=True)
    train_dataset = EnhancerDataset(
        args.train_pos,
        args.train_neg,
        args.sequence_length,
        augment=True,
        rc_prob=args.rc_prob,
        shift_prob=args.shift_prob,
        max_shift=args.max_shift,
    )
    validation_dataset = EnhancerDataset(
        args.validation_pos,
        args.validation_neg,
        args.sequence_length,
    )
    print(
        f"train: {len(train_dataset):,} records; "
        f"validation: {len(validation_dataset):,} records",
        flush=True,
    )

    generator = torch.Generator().manual_seed(args.seed)
    train_loader = build_loader(
        train_dataset, args.batch_size, args.num_workers, True, device, generator
    )
    validation_loader = build_loader(
        validation_dataset, args.batch_size, args.num_workers, False, device, generator
    )

    print(f"Loading and freezing AlphaGenome trunk from {args.weights}", flush=True)
    model = load_trunk(AlphaGenome(), str(args.weights), exclude_heads=True)
    for parameter in model.parameters():
        parameter.requires_grad = False
    model = remove_all_heads(model).to(device)
    model.eval()
    if any(parameter.requires_grad for parameter in model.parameters()):
        raise RuntimeError("AlphaGenome backbone contains trainable parameters")

    smoke_indices = (0, train_dataset.n_positive)
    smoke_sequences = torch.stack([train_dataset[index][0] for index in smoke_indices]).to(device)
    smoke_targets = torch.stack([train_dataset[index][1] for index in smoke_indices]).to(device)
    smoke_encoder = encode(model, smoke_sequences, device, use_amp)
    if smoke_encoder.shape[0] != 2 or smoke_encoder.shape[-1] != ENCODER_DIM:
        raise RuntimeError(f"Unexpected encoder shape: {tuple(smoke_encoder.shape)}")
    head = EnhancerHead(smoke_encoder.shape[1], args.hidden_size, args.dropout).to(device)
    smoke_logits = head(smoke_encoder.detach())
    smoke_loss = F.binary_cross_entropy_with_logits(smoke_logits.float(), smoke_targets.float())
    smoke_loss.backward()
    if smoke_logits.shape != (2,):
        raise RuntimeError(f"Unexpected logit shape: {tuple(smoke_logits.shape)}")
    if not any(parameter.grad is not None for parameter in head.parameters()):
        raise RuntimeError("EnhancerHead did not receive gradients")
    if any(parameter.grad is not None for parameter in model.parameters()):
        raise RuntimeError("Frozen AlphaGenome backbone received gradients")
    head.zero_grad(set_to_none=True)
    print(
        f"Preflight PASS: input={tuple(smoke_sequences.shape)}, "
        f"encoder={tuple(smoke_encoder.shape)}, logits={tuple(smoke_logits.shape)}",
        flush=True,
    )

    optimizer = torch.optim.Adam(
        head.parameters(), lr=args.learning_rate, weight_decay=args.weight_decay
    )
    scheduler = create_lr_scheduler(
        optimizer,
        warmup_steps=len(train_loader),
        total_steps=args.epochs * len(train_loader),
        schedule="cosine",
    )
    history: list[dict[str, float]] = []
    best_validation_auprc = -float("inf")
    patience_counter = 0

    for epoch in range(1, args.epochs + 1):
        train_metrics = train_epoch(
            model,
            head,
            train_loader,
            optimizer,
            scheduler,
            device,
            use_amp,
            args.threshold,
        )
        validation_metrics = evaluate(
            model,
            head,
            validation_loader,
            device,
            use_amp,
            args.threshold,
            "validation",
        )
        row: dict[str, float] = {"epoch": epoch}
        row.update({f"train_{key}": value for key, value in train_metrics.items()})
        row.update({f"validation_{key}": value for key, value in validation_metrics.items()})
        history.append(row)
        write_history(outputs["history"], history)

        improved = validation_metrics["auprc"] > best_validation_auprc
        print(
            f"epoch={epoch} train={train_metrics} validation={validation_metrics} "
            f"best={improved}",
            flush=True,
        )
        if improved:
            best_validation_auprc = validation_metrics["auprc"]
            patience_counter = 0
            save_checkpoint(
                outputs["checkpoint"],
                {
                    "epoch": epoch,
                    "head": head.state_dict(),
                    "validation_metrics": validation_metrics,
                    "config": config,
                    "encoder_positions": smoke_encoder.shape[1],
                    "encoder_dim": ENCODER_DIM,
                },
            )
        else:
            patience_counter += 1
            if patience_counter >= args.patience:
                print(f"Early stopping at epoch {epoch}", flush=True)
                break

    checkpoint = torch.load(outputs["checkpoint"], map_location=device, weights_only=True)
    head.load_state_dict(checkpoint["head"])
    restored_metrics = evaluate(
        model,
        head,
        validation_loader,
        device,
        use_amp,
        args.threshold,
        "validation",
    )
    validation_result = {
        "best_epoch": checkpoint["epoch"],
        "selection_metric": "validation_auprc",
        "validation": restored_metrics,
    }
    write_json(outputs["validation_metrics"], validation_result)

    print(
        "Loading held-out test FASTAs after validation-based model selection...",
        flush=True,
    )
    test_dataset = EnhancerDataset(
        args.test_pos,
        args.test_neg,
        args.sequence_length,
    )
    test_loader = build_loader(
        test_dataset, args.batch_size, args.num_workers, False, device, generator
    )
    test_metrics = evaluate(
        model,
        head,
        test_loader,
        device,
        use_amp,
        args.threshold,
        "test",
    )
    test_result = {
        "best_epoch": checkpoint["epoch"],
        "checkpoint_selection": "validation_auprc",
        "test": test_metrics,
        "elapsed_seconds": time.time() - started,
    }
    write_json(outputs["test_metrics"], test_result)
    print(json.dumps({**validation_result, **test_result}, indent=2, sort_keys=True), flush=True)
    print(f"Best head: {outputs['checkpoint']}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
