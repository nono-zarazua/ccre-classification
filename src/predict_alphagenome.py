#!/usr/bin/env python3
"""Score one positive/negative FASTA pair with a trained AlphaGenome head."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import torch
from torch.utils.data import DataLoader, Dataset

from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.transfer import load_trunk, remove_all_heads
from alphagenome_pytorch.extensions.finetuning.utils import sequence_to_onehot

from train_alphagenome import EnhancerHead, amp_context, normalize_encoder_output, read_fasta


class PredictionDataset(Dataset):
    def __init__(self, positive: Path, negative: Path, sequence_length: int) -> None:
        pos = read_fasta(positive)
        neg = read_fasta(negative)
        self.records = pos + neg
        self.labels = np.concatenate((np.ones(len(pos)), np.zeros(len(neg)))).astype(np.float32)
        self.sequence_length = sequence_length

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int):
        identifier, sequence = self.records[index]
        onehot = sequence_to_onehot(sequence).astype(np.float32)
        length = onehot.shape[0]
        if length < self.sequence_length:
            difference = self.sequence_length - length
            left = difference // 2
            onehot = np.pad(onehot, ((left, difference - left), (0, 0)))
        elif length > self.sequence_length:
            start = (length - self.sequence_length) // 2
            onehot = onehot[start : start + self.sequence_length]
        return identifier, torch.from_numpy(onehot), torch.tensor(self.labels[index])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("positive_fasta", type=Path)
    parser.add_argument("negative_fasta", type=Path)
    parser.add_argument("output_tsv", type=Path)
    parser.add_argument("--weights", type=Path, default=Path("weights/fold_0_weights.safetensors"))
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--num-workers", type=int, default=4)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for path in (args.checkpoint, args.positive_fasta, args.negative_fasta, args.weights):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")
    if args.output_tsv.exists():
        raise FileExistsError(f"Output already exists: {args.output_tsv}")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required for AlphaGenome inference")

    device = torch.device("cuda")
    checkpoint = torch.load(args.checkpoint, map_location="cpu", weights_only=False)
    config = checkpoint["config"]
    sequence_length = int(config["sequence_length"])
    dataset = PredictionDataset(args.positive_fasta, args.negative_fasta, sequence_length)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False,
                        num_workers=args.num_workers, pin_memory=True)

    trunk = load_trunk(AlphaGenome(), str(args.weights), exclude_heads=True)
    trunk = remove_all_heads(trunk).to(device).eval()
    for parameter in trunk.parameters():
        parameter.requires_grad = False

    head = EnhancerHead(int(checkpoint["encoder_positions"]),
                        int(config["hidden_size"]), float(config["dropout"]))
    head.load_state_dict(checkpoint["head"])
    head = head.to(device).eval()
    use_amp = not bool(config.get("no_amp", False))

    rows: list[tuple[str, int, float, float]] = []
    with torch.no_grad():
        for identifiers, sequences, labels in loader:
            sequences = sequences.to(device, non_blocking=True)
            organism = torch.zeros(sequences.shape[0], dtype=torch.long, device=device)
            with amp_context(device, use_amp):
                encoded = trunk(sequences, organism, encoder_only=True)["encoder_output"]
                logits = head(normalize_encoder_output(encoded)).float()
            probabilities = torch.sigmoid(logits)
            rows.extend(zip(identifiers, labels.int().tolist(), logits.cpu().tolist(),
                            probabilities.cpu().tolist()))

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output_tsv) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("sequence_id", "label", "logit", "probability"))
        writer.writerows(rows)
    temporary.replace(args.output_tsv)
    print(json.dumps({"output": str(args.output_tsv), "records": len(rows),
                      "positive": dataset.labels.sum().item(),
                      "negative": (dataset.labels == 0).sum().item(),
                      "checkpoint": str(args.checkpoint)}, indent=2))


if __name__ == "__main__":
    main()
