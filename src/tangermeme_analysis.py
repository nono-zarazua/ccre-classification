#!/usr/bin/env python3
"""Interpret a final AlphaGenome classifier with Tangermeme.

The primary workflow follows Tangermeme's recursive-seqlet/TOMTOM vignette:
attribute every held-out sequence, call signed recursive seqlets, annotate and
count every seqlet against a clustered MEME database, and draw representative
attribution logos. A separate strongest-site panel, balanced across sequence
class and attribution direction, reproduces the supervisor-style causal check
with repeated ablation, matched non-candidate controls, and bootstrap summaries.
Saturation mutagenesis remains an independent attribution sanity check.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import average_precision_score, roc_auc_score
from scipy.stats import fisher_exact

import alphagenome_pytorch
import tangermeme
from alphagenome_pytorch import AlphaGenome
from alphagenome_pytorch.extensions.finetuning.transfer import load_trunk, remove_all_heads
from alphagenome_pytorch.extensions.finetuning.utils import sequence_to_onehot
from tangermeme.ablate import ablate_annotations
from tangermeme.annotate import annotate_seqlets, count_annotations
from tangermeme.deep_lift_shap import BilinearOp, deep_lift_shap
from tangermeme.ersatz import dinucleotide_shuffle
from tangermeme.io import read_meme
from tangermeme.predict import predict
from tangermeme.plot import plot_logo
from tangermeme.saliency import saliency
from tangermeme.saturation_mutagenesis import saturation_mutagenesis
from tangermeme.seqlet import recursive_seqlets, tfmodisco_seqlets
from tangermeme.utils import set_seed


ENCODER_DIM = 1536
VALID_TASKS = {"evn", "svn", "evs"}


class AttributionHead(nn.Module):
    """Parameter-compatible copy of the trained binary head with module ReLU."""

    def __init__(self, n_positions: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        self.norm = nn.LayerNorm(ENCODER_DIM)
        self.fc1 = nn.Linear(n_positions * ENCODER_DIM, hidden_size)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(dropout)
        self.fc2 = nn.Linear(hidden_size, 1)

    def forward(self, encoder_output: torch.Tensor) -> torch.Tensor:
        if encoder_output.ndim != 3:
            raise ValueError(f"Expected 3D encoder output, got {tuple(encoder_output.shape)}")
        if encoder_output.shape[-1] != ENCODER_DIM:
            if encoder_output.shape[1] != ENCODER_DIM:
                raise ValueError(f"Cannot locate encoder dimension in {tuple(encoder_output.shape)}")
            encoder_output = encoder_output.transpose(1, 2)
        x = self.norm(encoder_output).flatten(1)
        return self.fc2(self.dropout(self.relu(self.fc1(x)))).squeeze(-1)


class TangermemeAlphaGenomeWrapper(nn.Module):
    """Map Tangermeme channel-first DNA to one raw classification logit."""

    def __init__(self, backbone: nn.Module, head: nn.Module, sequence_length: int) -> None:
        super().__init__()
        self.backbone = backbone
        self.head = head
        self.sequence_length = int(sequence_length)

    def forward(self, X: torch.Tensor) -> torch.Tensor:
        expected = (4, self.sequence_length)
        if X.ndim != 3 or tuple(X.shape[1:]) != expected:
            raise ValueError(f"Expected (batch, 4, {self.sequence_length}), got {tuple(X.shape)}")
        organism = torch.zeros(X.shape[0], dtype=torch.long, device=X.device)
        encoded = self.backbone(
            X.permute(0, 2, 1).contiguous(), organism, encoder_only=True
        )["encoder_output"]
        return self.head(encoded).unsqueeze(1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=sorted(VALID_TASKS), required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--weights", type=Path, required=True)
    parser.add_argument("--positive-fasta", type=Path, required=True)
    parser.add_argument("--negative-fasta", type=Path, required=True)
    parser.add_argument("--motif-meme", type=Path, required=True)
    parser.add_argument(
        "--motif-metadata",
        type=Path,
        required=True,
        help="HOCOMOCO clustered motif metadata (cluster_list.tsv)",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-per-class", type=int, default=128)
    parser.add_argument(
        "--all-sequences",
        action="store_true",
        help="Analyze every record in both FASTA files instead of sampling each class.",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--method", choices=("auto", "deep_lift_shap", "saliency"), default="auto")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--n-shuffles", type=int, default=10)
    parser.add_argument("--deeplift-tolerance", type=float, default=1e-3)
    parser.add_argument("--audit-sequences", type=int, default=4)
    parser.add_argument("--audit-shuffles", type=int, default=5)
    parser.add_argument("--tomtom-pvalue", type=float, default=0.05)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--candidate-window-size", type=int, default=21)
    parser.add_argument(
        "--seqlet-selection",
        choices=("fixed", "recursive", "modisco"),
        default="recursive",
        help=(
            "How attribution regions are called: the original fixed-width ranking, "
            "Tangermeme recursive seqlets, or the TF-MoDISco seqlet caller."
        ),
    )
    parser.add_argument("--candidate-fraction", type=float, default=0.02)
    parser.add_argument("--max-candidates-per-sign", type=int, default=3)
    parser.add_argument(
        "--max-candidate-sites",
        type=int,
        default=30000,
        help=(
            "Strongest class/sign-balanced sites reserved for ablation. All "
            "called seqlets are annotated and counted before this cap is applied; "
            "0 ablates every called seqlet."
        ),
    )
    parser.add_argument("--recursive-threshold", type=float, default=0.01)
    parser.add_argument("--recursive-min-length", type=int, default=4)
    parser.add_argument("--recursive-max-length", type=int, default=25)
    parser.add_argument("--recursive-flank", type=int, default=0)
    parser.add_argument("--modisco-window-size", type=int, default=21)
    parser.add_argument("--modisco-flank", type=int, default=10)
    parser.add_argument("--modisco-target-fdr", type=float, default=0.2)
    parser.add_argument("--modisco-min-passing-fraction", type=float, default=0.03)
    parser.add_argument("--modisco-max-passing-fraction", type=float, default=0.2)
    parser.add_argument(
        "--ablate-annotated-only",
        action="store_true",
        help="Ablate only candidates with a motif match (not recommended for primary runs).",
    )
    parser.add_argument("--controls-per-candidate", type=int, default=4)
    parser.add_argument("--min-known-fraction", type=float, default=0.8)
    parser.add_argument("--ablation-shuffles", type=int, default=10)
    parser.add_argument("--ablation-site-batch", type=int, default=64)
    parser.add_argument("--region-fdr", type=float, default=0.05)
    parser.add_argument("--control-quantile", type=float, default=0.95)
    parser.add_argument("--motif-fdr", type=float, default=0.05)
    parser.add_argument("--bootstrap-replicates", type=int, default=1000)
    parser.add_argument("--n-ism-sequences", type=int, default=8)
    parser.add_argument("--ism-batch-size", type=int, default=256)
    parser.add_argument(
        "--n-logo-sequences",
        type=int,
        default=8,
        help="Number of representative ablation-panel sequences to plot.",
    )
    parser.add_argument(
        "--logo-flank",
        type=int,
        default=50,
        help="Bases shown on either side of the representative seqlet.",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def validate_args(args: argparse.Namespace) -> None:
    for path in (
        args.checkpoint,
        args.weights,
        args.positive_fasta,
        args.negative_fasta,
        args.motif_meme,
        args.motif_metadata,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")
    for name in (
        "n_per_class", "batch_size", "n_shuffles", "audit_sequences",
        "audit_shuffles", "n_jobs", "candidate_window_size",
        "max_candidates_per_sign",
        "recursive_min_length", "recursive_max_length", "modisco_window_size",
        "controls_per_candidate", "ablation_shuffles", "ablation_site_batch",
        "bootstrap_replicates", "n_ism_sequences",
        "ism_batch_size", "n_logo_sequences",
    ):
        if getattr(args, name) <= 0:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if args.max_candidate_sites < 0:
        raise ValueError("--max-candidate-sites must be non-negative (0 disables the cap)")
    if args.recursive_flank < 0 or args.modisco_flank < 0:
        raise ValueError("Seqlet flanks must be non-negative")
    if args.logo_flank < 0:
        raise ValueError("--logo-flank must be non-negative")
    if args.recursive_min_length > args.recursive_max_length:
        raise ValueError("--recursive-min-length cannot exceed --recursive-max-length")
    if not 0 < args.recursive_threshold <= 1:
        raise ValueError("--recursive-threshold must be in (0, 1]")
    for name in (
        "candidate_fraction", "min_known_fraction", "region_fdr",
        "control_quantile", "motif_fdr", "modisco_target_fdr",
        "modisco_min_passing_fraction", "modisco_max_passing_fraction",
    ):
        if not 0 < getattr(args, name) <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be in (0, 1]")
    if not 0 < args.tomtom_pvalue <= 1:
        raise ValueError("--tomtom-pvalue must be in (0, 1]")
    if args.deeplift_tolerance <= 0:
        raise ValueError("--deeplift-tolerance must be positive")
    if args.modisco_min_passing_fraction > args.modisco_max_passing_fraction:
        raise ValueError(
            "--modisco-min-passing-fraction cannot exceed "
            "--modisco-max-passing-fraction"
        )
    if args.output_dir.exists() and not args.overwrite:
        raise FileExistsError(f"Output directory exists: {args.output_dir}")
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required; submit this analysis to a GPU node")


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
    if any(not sequence for _, sequence in records):
        raise ValueError(f"Empty sequence in {path}")
    identifiers = [identifier for identifier, _ in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate FASTA identifiers in {path}")
    return records


def center_adjust(onehot: np.ndarray, length: int) -> np.ndarray:
    observed = onehot.shape[0]
    if observed < length:
        difference = length - observed
        left = difference // 2
        return np.pad(onehot, ((left, difference - left), (0, 0)))
    if observed > length:
        start = (observed - length) // 2
        return onehot[start : start + length]
    return onehot


def sample_and_encode(
    positive_fasta: Path,
    negative_fasta: Path,
    n_per_class: int,
    sequence_length: int,
    seed: int,
    all_sequences: bool = False,
) -> tuple[torch.Tensor, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    classes = ((1, "positive", positive_fasta), (0, "negative", negative_fasta))
    tensors: list[torch.Tensor] = []
    rows: list[dict[str, Any]] = []
    for label, class_name, path in classes:
        records = read_fasta(path)
        if not all_sequences and n_per_class > len(records):
            raise ValueError(
                f"Requested {n_per_class:,} {class_name} records but {path} has {len(records):,}"
            )
        indices = (
            np.arange(len(records), dtype=int)
            if all_sequences
            else np.sort(rng.choice(len(records), size=n_per_class, replace=False))
        )
        for source_index in indices:
            identifier, sequence = records[int(source_index)]
            onehot = center_adjust(sequence_to_onehot(sequence).astype(np.float32), sequence_length)
            tensors.append(torch.from_numpy(onehot.T.astype(np.int8, copy=False)))
            rows.append(
                {
                    "example_idx": len(rows),
                    "sequence_id": f"{class_name}::{identifier}",
                    "fasta_id": identifier,
                    "label": label,
                    "class_name": class_name,
                    "source_fasta": str(path),
                    "source_index": int(source_index),
                    "original_length": len(sequence),
                    "analyzed_length": sequence_length,
                    "n_unknown_positions": int(np.all(onehot == 0, axis=1).sum()),
                }
            )
    return torch.stack(tensors), pd.DataFrame(rows)


def load_model(
    checkpoint_path: Path, weights_path: Path, device: torch.device
) -> tuple[TangermemeAlphaGenomeWrapper, dict[str, Any]]:
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    config = checkpoint["config"]
    sequence_length = int(config["sequence_length"])
    n_positions = int(checkpoint["encoder_positions"])
    backbone = remove_all_heads(
        load_trunk(AlphaGenome(), str(weights_path), exclude_heads=True)
    )
    for parameter in backbone.parameters():
        parameter.requires_grad_(False)
    head = AttributionHead(n_positions, int(config["hidden_size"]), float(config["dropout"]))
    missing, unexpected = head.load_state_dict(checkpoint["head"], strict=False)
    if missing or unexpected:
        raise RuntimeError(f"Head state mismatch; missing={missing}, unexpected={unexpected}")
    model = TangermemeAlphaGenomeWrapper(backbone, head, sequence_length).eval().to(device)
    return model, checkpoint


def native_logits(
    model: TangermemeAlphaGenomeWrapper,
    X: torch.Tensor,
    batch_size: int,
) -> torch.Tensor:
    """Predict through the native AlphaGenome path using explicit batches."""
    device = next(model.parameters()).device
    outputs = []
    with torch.no_grad():
        for start in range(0, len(X), batch_size):
            batch = X[start : start + batch_size].to(device).float()
            organism = torch.zeros(batch.shape[0], dtype=torch.long, device=device)
            encoded = model.backbone(
                batch.permute(0, 2, 1).contiguous(), organism, encoder_only=True
            )["encoder_output"]
            outputs.append(model.head(encoded).unsqueeze(1).cpu())
    return torch.cat(outputs)


def prediction_parity(
    model: TangermemeAlphaGenomeWrapper,
    X: torch.Tensor,
    device: torch.device,
    batch_size: int,
) -> tuple[torch.Tensor, dict[str, float]]:
    n = min(len(X), max(2, batch_size))
    native = native_logits(model, X[:n], batch_size)
    wrapped = predict(model, X[:n], batch_size=batch_size, device=device).cpu()
    difference = (native - wrapped).abs()
    result = {
        "n_sequences": n,
        "max_abs_difference": float(difference.max()),
        "mean_abs_difference": float(difference.mean()),
    }
    if result["max_abs_difference"] > 1e-5:
        raise RuntimeError(f"Native/wrapper prediction parity failed: {result}")
    logits = predict(model, X, batch_size=batch_size, device=device).cpu().squeeze(1)
    return logits, result


def deep_lift_rule_types() -> set[type[nn.Module]]:
    names = (
        "ReLU", "ReLU6", "RReLU", "SELU", "CELU", "GELU", "SiLU", "Mish",
        "GLU", "ELU", "LeakyReLU", "Sigmoid", "Tanh", "Softplus",
        "Softshrink", "LogSigmoid", "PReLU", "MaxPool1d", "MaxPool2d",
        "Softmax", "LayerNorm", "RMSNorm",
    )
    return {getattr(nn, name) for name in names if hasattr(nn, name)} | {BilinearOp}


LINEAR_OPS = {
    "linear", "conv1d", "conv2d", "conv3d", "conv_transpose1d", "add", "sub",
    "cat", "stack", "reshape", "view", "permute", "transpose", "flatten",
    "unsqueeze", "squeeze", "chunk", "split", "sum", "mean", "pad", "to",
    "type", "clone", "detach", "contiguous", "dim", "size", "dropout",
    "avg_pool1d", "avg_pool2d", "embedding", "batch_norm", "expand", "repeat",
    "narrow", "index_select", "roll", "__get__", "__getitem__",
    "_set_grad_enabled", "empty", "zeros", "ones", "arange", "t", "float",
    "half", "double", "cpu", "cuda", "numel", "shape", "is_floating_point",
    "masked_fill", "__setitem__", "neg", "flip",
}


def audit_operations(model: nn.Module, X: torch.Tensor) -> list[dict[str, str]]:
    ruled = deep_lift_rule_types()
    stack: list[nn.Module] = []
    found: set[tuple[str, str]] = set()

    def push(module: nn.Module, inputs: tuple[Any, ...]) -> None:
        del inputs
        stack.append(module)

    def pop(module: nn.Module, inputs: tuple[Any, ...], output: Any) -> None:
        del module, inputs, output
        stack.pop()

    class Trace(torch.overrides.TorchFunctionMode):
        def __torch_function__(self, func, types, args=(), kwargs=None):
            del types
            name = getattr(func, "__name__", "?")
            if name not in LINEAR_OPS and stack and type(stack[-1]) not in ruled:
                found.add((type(stack[-1]).__name__, name))
            return func(*args, **(kwargs or {}))

    hooks = [
        nn.modules.module.register_module_forward_pre_hook(push),
        nn.modules.module.register_module_forward_hook(pop),
    ]
    try:
        with Trace(), torch.no_grad():
            model(X.to(next(model.parameters()).device).float())
    finally:
        for hook in hooks:
            hook.remove()
    return [{"module": module, "operation": operation} for module, operation in sorted(found)]


def convergence_metrics(
    model: nn.Module,
    X: torch.Tensor,
    attributions: torch.Tensor,
    references: torch.Tensor,
    device: torch.device,
) -> dict[str, float]:
    model.eval()
    with torch.no_grad():
        y = model(X.to(device).float())[:, 0].float()
        ref_predictions = []
        for index in range(references.shape[1]):
            ref_predictions.append(model(references[:, index].to(device).float())[:, 0].float())
        y_reference = torch.stack(ref_predictions, dim=1).mean(dim=1)
    delta = attributions.sum(dim=(1, 2)).to(device) - (y - y_reference)
    scale = float((y - y_reference).abs().mean())
    mean_delta = float(delta.abs().mean())
    return {
        "mean_abs_delta": mean_delta,
        "max_abs_delta": float(delta.abs().max()),
        "output_difference_scale": scale,
        "relative_error": mean_delta / max(scale, 1e-12),
    }


def run_attribution_audit(
    model: nn.Module,
    X: torch.Tensor,
    args: argparse.Namespace,
    device: torch.device,
) -> tuple[str, dict[str, Any]]:
    audit_X = X[: min(args.audit_sequences, len(X))]
    payload: dict[str, Any] = {"requested_method": args.method}
    try:
        payload["unhooked_operations"] = audit_operations(model, audit_X[:1])
    except Exception as error:
        payload["operation_audit_error"] = f"{type(error).__name__}: {error}"
        payload["unhooked_operations"] = []

    references = dinucleotide_shuffle(
        audit_X,
        n=args.audit_shuffles,
        random_state=args.seed,
        allow_N=True,
        n_jobs=args.n_jobs,
    )
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            attr = deep_lift_shap(
                model,
                audit_X,
                target=0,
                batch_size=args.batch_size,
                references=references,
                n_shuffles=args.audit_shuffles,
                device=device,
                only_warn=True,
            )
        payload["deeplift"] = convergence_metrics(model, audit_X, attr, references, device)
        payload["deeplift"]["warnings"] = [str(item.message) for item in caught]
        payload["deeplift"]["passed"] = (
            payload["deeplift"]["relative_error"] <= args.deeplift_tolerance
        )
    except Exception as error:
        payload["deeplift"] = {
            "passed": False,
            "error": f"{type(error).__name__}: {error}",
        }

    if args.method == "deep_lift_shap" and not payload["deeplift"]["passed"]:
        raise RuntimeError(f"DeepLIFT audit failed: {payload['deeplift']}")
    selected = args.method
    if selected == "auto":
        selected = "deep_lift_shap" if payload["deeplift"]["passed"] else "saliency"
    payload["selected_method"] = selected
    payload["deeplift_tolerance"] = args.deeplift_tolerance
    return selected, payload


def attribute_sequences(
    model: nn.Module,
    X: torch.Tensor,
    method: str,
    args: argparse.Namespace,
    device: torch.device,
    seed: int,
) -> tuple[torch.Tensor, torch.Tensor | None, list[str]]:
    captured: list[str] = []
    if method == "deep_lift_shap":
        references = dinucleotide_shuffle(
            X,
            n=args.n_shuffles,
            random_state=seed,
            allow_N=True,
            n_jobs=args.n_jobs,
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            attr = deep_lift_shap(
                model,
                X,
                target=0,
                batch_size=args.batch_size,
                references=references,
                n_shuffles=args.n_shuffles,
                device=device,
                only_warn=True,
            )
        captured = [str(item.message) for item in caught]
        return attr.cpu(), references.cpu(), captured
    attr = saliency(
        model,
        X,
        target=0,
        batch_size=args.batch_size,
        gradient_correction=True,
        device=device,
    )
    return attr.cpu(), None, captured


def intervals_overlap(start: int, end: int, intervals: list[tuple[int, int]]) -> bool:
    return any(start < other_end and other_start < end for other_start, other_end in intervals)


def cap_candidates_stratified(
    candidates: pd.DataFrame,
    manifest: pd.DataFrame,
    maximum: int,
) -> pd.DataFrame:
    """Select the strongest candidates with equal class/sign quotas."""
    if maximum == 0 or len(candidates) <= maximum:
        return candidates.reset_index(drop=True)
    labels = manifest.set_index("example_idx")["label"]
    work = candidates.assign(
        _stratum_label=candidates["example_idx"].map(labels)
    )
    groups = list(work.groupby(["_stratum_label", "sign"], sort=True))
    quota = max(1, maximum // max(len(groups), 1))
    kept = []
    used: set[int] = set()
    for _, group in groups:
        chosen = group.nlargest(quota, "absolute_attribution")
        kept.append(chosen)
        used.update(chosen.index.tolist())
    selected = pd.concat(kept) if kept else work.iloc[0:0]
    remaining = maximum - len(selected)
    if remaining > 0:
        extra = work.loc[~work.index.isin(used)].nlargest(remaining, "absolute_attribution")
        selected = pd.concat((selected, extra))
    selected = selected.drop(columns="_stratum_label").nlargest(
        maximum, "absolute_attribution"
    )
    return selected.reset_index(drop=True)


def select_fixed_candidates(
    X: torch.Tensor,
    attribution_scores: torch.Tensor,
    manifest: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Select strong, signed, non-overlapping fixed-width saliency windows."""
    width = args.candidate_window_size
    if width > attribution_scores.shape[1]:
        raise ValueError("--candidate-window-size exceeds the analyzed sequence length")
    kernel = np.ones(width, dtype=np.float64)
    rows: list[dict[str, Any]] = []
    for example_idx in range(len(X)):
        score = attribution_scores[example_idx].numpy().astype(np.float64, copy=False)
        rolling = np.convolve(score, kernel, mode="valid")
        known = (X[example_idx].sum(dim=0) > 0).numpy().astype(np.float64)
        known_fraction = np.convolve(known, kernel, mode="valid") / width
        eligible = np.flatnonzero(known_fraction >= args.min_known_fraction)
        n_to_take = min(
            args.max_candidates_per_sign,
            max(1, int(np.ceil(len(eligible) * args.candidate_fraction))),
        )
        occupied: list[tuple[int, int]] = []
        for sign, direction in (("positive", 1.0), ("negative", -1.0)):
            ranked = eligible[np.argsort(direction * rolling[eligible])[::-1]]
            taken = 0
            for start_value in ranked:
                start = int(start_value)
                end = start + width
                signed_strength = float(direction * rolling[start])
                if signed_strength <= 0:
                    break
                if intervals_overlap(start, end, occupied):
                    continue
                attribution = float(rolling[start])
                rows.append(
                    {
                        "example_idx": example_idx,
                        "start": start,
                        "end": end,
                        "sign": sign,
                        "attribution": attribution,
                        "absolute_attribution": abs(attribution),
                        "known_fraction": float(known_fraction[start]),
                    }
                )
                occupied.append((start, end))
                taken += 1
                if taken >= n_to_take:
                    break
    candidates = pd.DataFrame(rows)
    if candidates.empty:
        raise RuntimeError("Saliency screening produced no eligible signed candidate windows")
    candidates.insert(0, "screen_candidate_id", np.arange(len(candidates), dtype=int))
    return candidates


def _known_fraction(X: torch.Tensor, example_idx: int, start: int, end: int) -> float:
    return float((X[example_idx, :, start:end].sum(dim=0) > 0).float().mean())


def _normalize_seqlet_candidates(
    X: torch.Tensor,
    seqlets: pd.DataFrame,
    min_known_fraction: float,
) -> pd.DataFrame:
    """Validate Tangermeme seqlets and add columns used downstream."""
    required = {"example_idx", "start", "end", "attribution"}
    missing = required.difference(seqlets.columns)
    if missing:
        raise ValueError(f"Seqlet caller omitted required columns: {sorted(missing)}")
    result = seqlets.loc[:, ["example_idx", "start", "end", "attribution"]].copy()
    result = result.astype(
        {"example_idx": "int64", "start": "int64", "end": "int64", "attribution": "float64"}
    )
    valid = (
        result["example_idx"].between(0, len(X) - 1)
        & (result["start"] >= 0)
        & (result["end"] > result["start"])
        & (result["end"] <= X.shape[-1])
        & (result["attribution"] != 0)
    )
    result = result.loc[valid].copy()
    result["sign"] = np.where(result["attribution"] > 0, "positive", "negative")
    result["absolute_attribution"] = result["attribution"].abs()
    result["known_fraction"] = [
        _known_fraction(X, int(row.example_idx), int(row.start), int(row.end))
        for row in result.itertuples(index=False)
    ]
    result = result.loc[result["known_fraction"] >= min_known_fraction].copy()
    return result.reset_index(drop=True)


def select_recursive_candidates(
    X: torch.Tensor,
    attribution_scores: torch.Tensor,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Call positive and negative recursive seqlets and resolve sign overlaps."""
    kwargs = {
        "threshold": args.recursive_threshold,
        "min_seqlet_len": args.recursive_min_length,
        "max_seqlet_len": args.recursive_max_length,
        "additional_flanks": args.recursive_flank,
    }
    positive = recursive_seqlets(attribution_scores, **kwargs).copy()
    negative = recursive_seqlets(-attribution_scores, **kwargs).copy()
    negative["attribution"] = -negative["attribution"]

    # Tangermeme recommends dropping negative calls overlapping positive calls.
    positive_intervals = {
        int(example_idx): list(zip(group["start"].astype(int), group["end"].astype(int)))
        for example_idx, group in positive.groupby("example_idx")
    }
    keep_negative = [
        not intervals_overlap(
            int(row.start), int(row.end),
            positive_intervals.get(int(row.example_idx), []),
        )
        for row in negative.itertuples(index=False)
    ]
    negative = negative.loc[keep_negative]
    candidates = pd.concat((positive, negative), ignore_index=True, sort=False)
    return _normalize_seqlet_candidates(X, candidates, args.min_known_fraction)


def select_modisco_candidates(
    X: torch.Tensor,
    attribution_scores: torch.Tensor,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Call signed regions with Tangermeme's TF-MoDISco seqlet procedure."""
    seqlets = tfmodisco_seqlets(
        attribution_scores,
        window_size=args.modisco_window_size,
        flank=args.modisco_flank,
        target_fdr=args.modisco_target_fdr,
        min_passing_frac=args.modisco_min_passing_fraction,
        max_passing_frac=args.modisco_max_passing_fraction,
    )
    return _normalize_seqlet_candidates(X, seqlets, args.min_known_fraction)


def select_candidates(
    X: torch.Tensor,
    attribution_scores: torch.Tensor,
    manifest: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    if args.seqlet_selection == "fixed":
        candidates = select_fixed_candidates(X, attribution_scores, manifest, args)
    elif args.seqlet_selection == "recursive":
        candidates = select_recursive_candidates(X, attribution_scores, args)
    else:
        candidates = select_modisco_candidates(X, attribution_scores, args)
    if candidates.empty:
        raise RuntimeError(f"{args.seqlet_selection} seqlet selection produced no candidates")
    if "screen_candidate_id" in candidates:
        candidates = candidates.drop(columns="screen_candidate_id")
    candidates.insert(0, "screen_candidate_id", np.arange(len(candidates), dtype=int))
    return candidates


def select_matched_controls(
    X: torch.Tensor,
    attribution_scores: torch.Tensor,
    candidates: pd.DataFrame,
    manifest: pd.DataFrame,
    args: argparse.Namespace,
    excluded_candidates: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Draw random non-candidate controls of equal length.

    Prefer the candidate's own sequence. If its known sequence is fully occupied
    by candidates, fall back to another held-out sequence from the same class.
    Controls are not selected for low saliency because conditioning the null on
    low saliency would make candidate/control separation partly circular.
    """
    rng = np.random.default_rng(args.seed + 300_000)
    labels = manifest.set_index("example_idx")["label"].to_dict()
    examples_by_label = {
        int(label): group["example_idx"].astype(int).tolist()
        for label, group in manifest.groupby("label", sort=True)
    }
    excluded = candidates if excluded_candidates is None else excluded_candidates
    candidates_by_example = {
        int(example_idx): list(zip(group["start"].astype(int), group["end"].astype(int)))
        for example_idx, group in excluded.groupby("example_idx")
    }
    used_control_intervals: dict[int, list[tuple[int, int]]] = {
        int(example_idx): [] for example_idx in manifest["example_idx"]
    }

    def available_starts(
        example_idx: int, width: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        kernel = np.ones(width, dtype=np.float64)
        score = attribution_scores[example_idx].numpy().astype(np.float64, copy=False)
        rolling = np.convolve(score, kernel, mode="valid")
        known = (X[example_idx].sum(dim=0) > 0).numpy().astype(np.float64)
        known_fraction = np.convolve(known, kernel, mode="valid") / width
        candidate_intervals = candidates_by_example.get(example_idx, [])
        available = np.asarray(
            [
                start for start in np.flatnonzero(known_fraction >= args.min_known_fraction)
                if not intervals_overlap(
                    int(start), int(start) + width, candidate_intervals
                )
                and not intervals_overlap(
                    int(start), int(start) + width, used_control_intervals[example_idx]
                )
            ],
            dtype=int,
        )
        return available, rolling, known_fraction

    rows: list[dict[str, Any]] = []
    for candidate in candidates.itertuples(index=False):
        candidate_example_idx = int(candidate.example_idx)
        width = int(candidate.end) - int(candidate.start)
        label = int(labels[candidate_example_idx])
        alternatives = [
            index for index in examples_by_label[label] if index != candidate_example_idx
        ]
        rng.shuffle(alternatives)
        search_order = [candidate_example_idx] + alternatives
        for _ in range(args.controls_per_candidate):
            selected: tuple[int, int, np.ndarray, np.ndarray] | None = None
            for control_example_idx in search_order:
                available, rolling, known_fraction = available_starts(
                    control_example_idx, width
                )
                if len(available):
                    selected = (
                        control_example_idx,
                        int(rng.choice(available)),
                        rolling,
                        known_fraction,
                    )
                    break
            if selected is None:
                raise RuntimeError(
                    "No class-matched non-candidate control window available for "
                    f"candidate {candidate.candidate_id}"
                )
            control_example_idx, start, rolling, known_fraction = selected
            end = start + width
            rows.append(
                {
                    "control_id": len(rows),
                    "matched_candidate_id": int(candidate.candidate_id),
                    "candidate_example_idx": candidate_example_idx,
                    "example_idx": control_example_idx,
                    "control_match": (
                        "same_sequence"
                        if control_example_idx == candidate_example_idx
                        else "same_class"
                    ),
                    "start": start,
                    "end": end,
                    "sign": candidate.sign,
                    "attribution": float(rolling[start]),
                    "absolute_attribution": abs(float(rolling[start])),
                    "known_fraction": float(known_fraction[start]),
                }
            )
            used_control_intervals[control_example_idx].append((start, end))
    return pd.DataFrame(rows)


def annotate_windows(
    X: torch.Tensor,
    windows: pd.DataFrame,
    motif_path: Path,
    motif_names: list[str],
    pvalue_threshold: float,
    n_jobs: int,
) -> pd.DataFrame:
    result = windows.copy()
    if result.empty:
        result["motif_index"] = pd.Series(dtype="int64")
        result["tomtom_pvalue"] = pd.Series(dtype="float64")
        result["motif"] = pd.Series(dtype="object")
        return result

    # ``annotate_seqlets`` uses the first three columns positionally as
    # (example_idx, start, end).  Keep analysis identifiers such as
    # candidate_id out of that positional interface.
    coordinates = result.loc[:, ["example_idx", "start", "end"]].copy()
    coordinates = coordinates.astype(
        {"example_idx": "int64", "start": "int64", "end": "int64"}
    )
    if not coordinates["example_idx"].between(0, len(X) - 1).all():
        invalid = coordinates.loc[
            ~coordinates["example_idx"].between(0, len(X) - 1), "example_idx"
        ]
        raise ValueError(
            "Motif-annotation example indices fall outside the sequence tensor: "
            f"min={invalid.min()}, max={invalid.max()}, n_examples={len(X)}"
        )
    if not (
        (coordinates["start"] >= 0)
        & (coordinates["end"] > coordinates["start"])
        & (coordinates["end"] <= X.shape[-1])
    ).all():
        raise ValueError("Motif-annotation coordinates fall outside sequence bounds")

    indices, pvalues = annotate_seqlets(
        X,
        coordinates,
        str(motif_path),
        n_nearest=1,
        n_jobs=n_jobs,
    )
    result["motif_index"] = indices[:, 0].cpu().numpy().astype(int)
    result["tomtom_pvalue"] = pvalues[:, 0].cpu().numpy()
    accepted = (result["motif_index"] >= 0) & (result["tomtom_pvalue"] <= pvalue_threshold)
    result["motif"] = pd.NA
    result.loc[accepted, "motif"] = [
        motif_names[index] for index in result.loc[accepted, "motif_index"]
    ]
    return result


def summarize_all_annotations(
    seqlets: pd.DataFrame,
    manifest: pd.DataFrame,
    motif_names: list[str],
    cluster_metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Count every accepted TOMTOM annotation before ablation selection.

    The sparse per-sequence table is convenient for downstream enrichment,
    while the motif/sign table gives a complete descriptive inventory.  The
    Tangermeme counter is used with int32 explicitly to avoid uint8 overflow.
    """
    accepted = seqlets.loc[
        seqlets["motif"].notna() & (seqlets["motif_index"] >= 0)
    ].copy()
    per_sequence_columns = [
        "example_idx", "sequence_id", "label", "class_name", "motif_index",
        "motif", "sign", "n_sites", "mean_attribution", "total_attribution",
    ]
    summary_columns = [
        "motif_index", "motif", "representative_motif", "representative_gene",
        "cluster_size", "primary_family", "primary_subfamily",
        "clustered_gene_symbols", "sign", "n_sites", "n_sequences",
        "all_signs_n_sites", "all_signs_n_sequences",
        "all_signs_median_attribution", "all_signs_attribution_q25",
        "all_signs_attribution_q75", "all_signs_attribution_q90",
        "all_signs_attribution_q95", "n_positive_sites", "n_negative_sites",
        "n_positive_sequences", "n_negative_sequences", "mean_attribution",
        "median_attribution", "attribution_q25", "attribution_q75",
        "attribution_q90", "attribution_q95", "min_attribution",
        "max_attribution", "total_attribution", "mean_tomtom_pvalue",
    ]
    if accepted.empty:
        return (
            pd.DataFrame(columns=per_sequence_columns),
            pd.DataFrame(columns=summary_columns),
        )

    accepted["positive_example_idx"] = accepted["example_idx"].where(
        accepted["label"] == 1
    )
    accepted["negative_example_idx"] = accepted["example_idx"].where(
        accepted["label"] == 0
    )
    annotation_pairs = torch.as_tensor(
        accepted[["example_idx", "motif_index"]].to_numpy(copy=True),
        dtype=torch.long,
    )
    all_sign_counts = count_annotations(
        annotation_pairs,
        dtype=torch.int32,
        shape=(len(manifest), len(motif_names)),
        dim=0,
    ).numpy()

    identity = manifest[["example_idx", "sequence_id", "label", "class_name"]]
    per_sequence = (
        accepted.groupby(
            ["example_idx", "motif_index", "motif", "sign"],
            observed=True,
            dropna=False,
        )
        .agg(
            n_sites=("motif", "size"),
            mean_attribution=("attribution", "mean"),
            total_attribution=("attribution", "sum"),
        )
        .reset_index()
        .merge(identity, on="example_idx", how="left")
    )
    per_sequence = per_sequence[per_sequence_columns]

    summary = (
        accepted.groupby(["motif_index", "motif", "sign"], observed=True)
        .agg(
            n_sites=("motif", "size"),
            n_sequences=("example_idx", "nunique"),
            n_positive_sites=("label", lambda values: int((values == 1).sum())),
            n_negative_sites=("label", lambda values: int((values == 0).sum())),
            n_positive_sequences=("positive_example_idx", "nunique"),
            n_negative_sequences=("negative_example_idx", "nunique"),
            mean_attribution=("attribution", "mean"),
            median_attribution=("attribution", "median"),
            attribution_q25=("attribution", lambda values: values.quantile(0.25)),
            attribution_q75=("attribution", lambda values: values.quantile(0.75)),
            attribution_q90=("attribution", lambda values: values.quantile(0.90)),
            attribution_q95=("attribution", lambda values: values.quantile(0.95)),
            min_attribution=("attribution", "min"),
            max_attribution=("attribution", "max"),
            total_attribution=("attribution", "sum"),
            mean_tomtom_pvalue=("tomtom_pvalue", "mean"),
        )
        .reset_index()
    )
    summary["all_signs_n_sites"] = summary["motif_index"].map(
        lambda index: int(all_sign_counts[int(index)])
    )
    all_signs = accepted.groupby("motif_index")["attribution"].agg(
        all_signs_n_sequences=lambda values: accepted.loc[
            values.index, "example_idx"
        ].nunique(),
        all_signs_median_attribution="median",
        all_signs_attribution_q25=lambda values: values.quantile(0.25),
        all_signs_attribution_q75=lambda values: values.quantile(0.75),
        all_signs_attribution_q90=lambda values: values.quantile(0.90),
        all_signs_attribution_q95=lambda values: values.quantile(0.95),
    )
    for column in all_signs.columns:
        summary[column] = summary["motif_index"].map(all_signs[column])
    summary["representative_motif"] = summary["motif"].str.split().str[0]
    metadata = cluster_metadata.rename(
        columns={
            "Representative_Motif": "representative_motif",
            "Representative_Motif_Gene_Symbol": "representative_gene",
            "Cluster_Size": "cluster_size",
            "Primary_Family": "primary_family",
            "Primary_Subfamily": "primary_subfamily",
            "Clustered_Motifs_Gene_Symbols": "clustered_gene_symbols",
        }
    )
    summary = summary.merge(
        metadata[
            [
                "representative_motif", "representative_gene", "cluster_size",
                "primary_family", "primary_subfamily", "clustered_gene_symbols",
            ]
        ],
        on="representative_motif",
        how="left",
        validate="many_to_one",
    )
    if summary["cluster_size"].isna().any():
        missing = sorted(summary.loc[summary["cluster_size"].isna(), "representative_motif"].unique())
        raise ValueError(f"Accepted motifs absent from cluster metadata: {missing}")
    summary = summary[summary_columns].sort_values(
        ["all_signs_n_sites", "n_sites", "motif", "sign"],
        ascending=[False, False, True, True],
    )
    return per_sequence, summary


def select_representative_logo_sites(
    candidates: pd.DataFrame,
    maximum: int,
) -> pd.DataFrame:
    """Choose unique, annotated examples while balancing class and sign."""
    if maximum <= 0:
        return candidates.iloc[0:0].copy()
    annotated = candidates.loc[candidates["motif"].notna()].copy()
    if annotated.empty:
        return annotated
    groups = {
        key: group.sort_values("absolute_attribution", ascending=False)
        for key, group in annotated.groupby(["label", "sign"], sort=True)
    }
    positions = {key: 0 for key in groups}
    selected = []
    used_examples: set[int] = set()
    while len(selected) < maximum:
        progress = False
        for key in sorted(groups):
            group = groups[key]
            position = positions[key]
            while position < len(group):
                row = group.iloc[position]
                position += 1
                if int(row["example_idx"]) in used_examples:
                    continue
                selected.append(row)
                used_examples.add(int(row["example_idx"]))
                progress = True
                break
            positions[key] = position
            if len(selected) >= maximum:
                break
        if not progress:
            break
    return pd.DataFrame(selected).reset_index(drop=True)


def write_attribution_logos(
    attribution: torch.Tensor,
    all_seqlets: pd.DataFrame,
    ablation_candidates: pd.DataFrame,
    manifest: pd.DataFrame,
    args: argparse.Namespace,
    output_dir: Path,
) -> pd.DataFrame:
    """Render representative attribution logos with TOMTOM motif labels."""
    logo_sites = select_representative_logo_sites(
        ablation_candidates, args.n_logo_sequences
    )
    columns = [
        "logo_rank", "example_idx", "sequence_id", "label", "class_name",
        "sign", "motif", "start", "end", "absolute_attribution",
        "plot_start", "plot_end", "png", "svg",
    ]
    if logo_sites.empty:
        return pd.DataFrame(columns=columns)

    logo_dir = output_dir / "attribution_logos"
    logo_dir.mkdir(parents=True, exist_ok=True)
    identity = manifest.set_index("example_idx")
    rows = []
    sequence_length = attribution.shape[-1]
    for rank, site in enumerate(logo_sites.itertuples(index=False), start=1):
        example_idx = int(site.example_idx)
        plot_start = max(0, int(site.start) - args.logo_flank)
        plot_end = min(sequence_length, int(site.end) + args.logo_flank)
        annotations = all_seqlets.loc[
            (all_seqlets["example_idx"] == example_idx)
            & all_seqlets["motif"].notna()
            & (all_seqlets["start"] > plot_start)
            & (all_seqlets["end"] < plot_end)
        ].copy()
        plot_annotations = annotations[
            ["motif", "start", "end", "attribution"]
        ].rename(columns={"motif": "motif_name"})

        fig, ax = plt.subplots(figsize=(13, 3.5))
        plot_logo(
            attribution[example_idx],
            ax=ax,
            start=plot_start,
            end=plot_end,
            annotations=plot_annotations if not plot_annotations.empty else None,
            score_key="attribution",
            show_score=False,
            n_tracks=4,
        )
        record = identity.loc[example_idx]
        ax.set_title(
            f"{args.task.upper()} {record['class_name']} | {site.sign} attribution | "
            f"{site.motif}"
        )
        ax.set_xlabel("Position in analyzed sequence")
        ax.set_ylabel("Attribution")
        fig.tight_layout()
        stem = f"logo_{rank:02d}_{record['class_name']}_{site.sign}_example_{example_idx}"
        png = logo_dir / f"{stem}.png"
        svg = logo_dir / f"{stem}.svg"
        fig.savefig(png, dpi=200, bbox_inches="tight")
        fig.savefig(svg, bbox_inches="tight")
        plt.close(fig)
        rows.append(
            {
                "logo_rank": rank,
                "example_idx": example_idx,
                "sequence_id": record["sequence_id"],
                "label": int(record["label"]),
                "class_name": record["class_name"],
                "sign": site.sign,
                "motif": site.motif,
                "start": int(site.start),
                "end": int(site.end),
                "absolute_attribution": float(site.absolute_attribution),
                "plot_start": plot_start,
                "plot_end": plot_end,
                "png": str(png.relative_to(output_dir)),
                "svg": str(svg.relative_to(output_dir)),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def run_ablation(
    model: nn.Module,
    X: torch.Tensor,
    sites: pd.DataFrame,
    args: argparse.Namespace,
    device: torch.device,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_columns = list(sites.columns) + [
        "native_logit", "mean_ablated_logit", "mean_ablation_delta",
        "median_ablation_delta", "sd_ablation_delta",
    ]
    replicate_columns = ["site_id", "ablation_replicate", "ablated_logit", "ablation_delta"]
    if sites.empty:
        return pd.DataFrame(columns=summary_columns), pd.DataFrame(columns=replicate_columns)
    summaries = []
    replicate_rows = []
    for chunk_start in range(0, len(sites), args.ablation_site_batch):
        chunk = sites.iloc[chunk_start : chunk_start + args.ablation_site_batch].copy()
        annotations = torch.tensor(
            chunk[["example_idx", "start", "end"]].to_numpy(), dtype=torch.long
        )
        result = ablate_annotations(
            model,
            X,
            annotations,
            n=args.ablation_shuffles,
            random_state=args.seed + chunk_start,
            batch_size=args.batch_size,
            device=device,
        )
        n_sites = len(chunk)
        before = result.y_befores.reshape(n_sites, -1).mean(dim=1).cpu().numpy()
        after = result.y_afters.reshape(n_sites, args.ablation_shuffles, -1).mean(dim=2).cpu().numpy()
        deltas = before[:, None] - after
        chunk["native_logit"] = before
        chunk["mean_ablated_logit"] = after.mean(axis=1)
        chunk["mean_ablation_delta"] = deltas.mean(axis=1)
        chunk["median_ablation_delta"] = np.median(deltas, axis=1)
        chunk["sd_ablation_delta"] = deltas.std(axis=1, ddof=1) if args.ablation_shuffles > 1 else 0.0
        summaries.append(chunk)
        for local_idx, site_id in enumerate(chunk["site_id"].astype(str)):
            for replicate in range(args.ablation_shuffles):
                replicate_rows.append(
                    {
                        "site_id": site_id,
                        "ablation_replicate": replicate,
                        "ablated_logit": float(after[local_idx, replicate]),
                        "ablation_delta": float(deltas[local_idx, replicate]),
                    }
                )
    return pd.concat(summaries, ignore_index=True), pd.DataFrame(replicate_rows)


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    finite = np.flatnonzero(np.isfinite(values))
    if len(finite) == 0:
        return adjusted
    order = finite[np.argsort(values[finite])]
    ranked = values[order]
    corrected = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    corrected = np.minimum.accumulate(corrected[::-1])[::-1]
    adjusted[order] = np.minimum(corrected, 1.0)
    return adjusted


def test_candidate_regions(
    candidate_ablation: pd.DataFrame,
    control_ablation: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Test candidate magnitudes against pooled class/sign-matched controls."""
    null_rows = []
    lookup: dict[tuple[int, str], np.ndarray] = {}
    for (label, sign), group in control_ablation.groupby(["label", "sign"], sort=True):
        null = group["mean_ablation_delta"].abs().to_numpy(dtype=float)
        lookup[(int(label), str(sign))] = null
        threshold = float(np.quantile(null, args.control_quantile))
        for value in null:
            null_rows.append(
                {
                    "label": int(label),
                    "sign": str(sign),
                    "absolute_control_effect": float(value),
                    "control_quantile_threshold": threshold,
                }
            )
    tested = candidate_ablation.copy()
    pvalues = []
    thresholds = []
    null_sizes = []
    for row in tested.itertuples(index=False):
        null = lookup.get((int(row.label), str(row.sign)))
        if null is None or len(null) == 0:
            raise RuntimeError(f"No matched null for label={row.label}, sign={row.sign}")
        effect = abs(float(row.mean_ablation_delta))
        pvalues.append((1 + int(np.sum(null >= effect))) / (len(null) + 1))
        thresholds.append(float(np.quantile(null, args.control_quantile)))
        null_sizes.append(len(null))
    tested["empirical_pvalue"] = pvalues
    tested["adjusted_pvalue"] = benjamini_hochberg(np.asarray(pvalues))
    tested["control_effect_threshold"] = thresholds
    tested["matched_null_size"] = null_sizes
    tested["important_region"] = (
        (tested["adjusted_pvalue"] <= args.region_fdr)
        & (tested["mean_ablation_delta"].abs() > tested["control_effect_threshold"])
    )
    return tested, pd.DataFrame(null_rows)


def read_cluster_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t")
    required = {
        "Representative_Motif", "Cluster_Size", "Clustered_Motifs",
        "Primary_Family", "Primary_Subfamily",
        "Representative_Motif_Gene_Symbol", "Clustered_Motifs_Gene_Symbols",
    }
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Motif metadata is missing columns: {sorted(missing)}")
    if metadata["Representative_Motif"].duplicated().any():
        raise ValueError("Motif metadata contains duplicate representative motifs")
    return metadata


def summarize_clustered_motifs(
    tested: pd.DataFrame,
    manifest: pd.DataFrame,
    cluster_metadata: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Summarize all annotated top-ranked sites without post-selection testing.

    Enrichment is tested at sequence level between labels within each saliency
    sign. Ablation effects are descriptive causal estimates with bootstrap
    intervals; they are not retested after selecting important regions.
    """
    annotated = tested.loc[tested["motif"].notna()].copy()
    columns = [
        "motif", "representative_motif", "representative_gene", "cluster_size",
        "primary_family", "primary_subfamily", "clustered_gene_symbols",
        "sign", "n_sites", "n_sequences", "n_important_regions", "mean_attribution",
        "total_attribution", "mean_tomtom_pvalue", "mean_ablation_delta",
        "median_ablation_delta", "sign_consistency", "bootstrap_ci_low",
        "bootstrap_ci_high", "enrichment_odds_ratio", "enrichment_pvalue",
        "enrichment_adjusted_pvalue",
        "n_positive_sequences", "n_negative_sequences",
        "positive_sequence_fraction", "negative_sequence_fraction",
        "enriched_cluster",
    ]
    if annotated.empty:
        return pd.DataFrame(columns=columns)
    rng = np.random.default_rng(args.seed + 400_000)
    metadata = cluster_metadata.set_index("Representative_Motif")
    rows = []
    candidate_sequences = {
        (str(sign), int(label)): group["example_idx"].nunique()
        for (sign, label), group in tested.groupby(["sign", "label"], sort=True)
    }
    for (motif, sign), group in annotated.groupby(["motif", "sign"], sort=True):
        effects = group["mean_ablation_delta"].to_numpy(dtype=float)
        bootstrap = np.asarray([
            rng.choice(effects, size=len(effects), replace=True).mean()
            for _ in range(args.bootstrap_replicates)
        ])
        positive_sequences = group.loc[group["label"] == 1, "example_idx"].nunique()
        negative_sequences = group.loc[group["label"] == 0, "example_idx"].nunique()
        n_positive = candidate_sequences.get((str(sign), 1), 0)
        n_negative = candidate_sequences.get((str(sign), 0), 0)
        table = [
            [positive_sequences, n_positive - positive_sequences],
            [negative_sequences, n_negative - negative_sequences],
        ]
        odds_ratio, enrichment_pvalue = fisher_exact(table, alternative="two-sided")
        motif_id = str(motif).split()[0]
        if motif_id not in metadata.index:
            raise ValueError(f"Motif is absent from cluster metadata: {motif_id}")
        cluster = metadata.loc[motif_id]
        rows.append(
            {
                "motif": motif,
                "representative_motif": motif_id,
                "representative_gene": cluster["Representative_Motif_Gene_Symbol"],
                "cluster_size": int(cluster["Cluster_Size"]),
                "primary_family": cluster["Primary_Family"],
                "primary_subfamily": cluster["Primary_Subfamily"],
                "clustered_gene_symbols": cluster["Clustered_Motifs_Gene_Symbols"],
                "sign": sign,
                "n_sites": len(group),
                "n_sequences": group["example_idx"].nunique(),
                "n_important_regions": int(group["important_region"].sum()),
                "mean_attribution": group["attribution"].mean(),
                "total_attribution": group["attribution"].sum(),
                "mean_tomtom_pvalue": group["tomtom_pvalue"].mean(),
                "mean_ablation_delta": effects.mean(),
                "median_ablation_delta": np.median(effects),
                "sign_consistency": max(np.mean(effects > 0), np.mean(effects < 0)),
                "bootstrap_ci_low": np.quantile(bootstrap, 0.025),
                "bootstrap_ci_high": np.quantile(bootstrap, 0.975),
                "enrichment_odds_ratio": odds_ratio,
                "enrichment_pvalue": enrichment_pvalue,
                "n_positive_sequences": positive_sequences,
                "n_negative_sequences": negative_sequences,
                "positive_sequence_fraction": positive_sequences / n_positive if n_positive else np.nan,
                "negative_sequence_fraction": negative_sequences / n_negative if n_negative else np.nan,
            }
        )
    summary = pd.DataFrame(rows)
    summary["enrichment_adjusted_pvalue"] = benjamini_hochberg(
        summary["enrichment_pvalue"].to_numpy()
    )
    summary["enriched_cluster"] = summary["enrichment_adjusted_pvalue"] <= args.motif_fdr
    return summary[columns].sort_values(
        ["enriched_cluster", "enrichment_adjusted_pvalue", "n_sites"],
        ascending=[False, True, False],
    )


def run_ism_validation(
    model: nn.Module,
    X: torch.Tensor,
    attribution_scores: torch.Tensor,
    candidates: pd.DataFrame,
    manifest: pd.DataFrame,
    args: argparse.Namespace,
    device: torch.device,
) -> tuple[torch.Tensor, pd.DataFrame, list[int]]:
    importance = attribution_scores.abs().max(dim=1).values
    candidate_rank = []
    if not candidates.empty:
        candidate_rank = (
            candidates.assign(abs_attr=candidates["attribution"].abs())
            .groupby("example_idx")["abs_attr"]
            .max()
            .sort_values(ascending=False)
            .index.tolist()
        )
    global_rank = list(dict.fromkeys(
        [int(index) for index in candidate_rank]
        + torch.argsort(importance, descending=True).tolist()
    ))
    positive_target = (args.n_ism_sequences + 1) // 2
    negative_target = args.n_ism_sequences // 2
    labels = manifest.set_index("example_idx")["label"].to_dict()
    positive = [index for index in global_rank if labels[index] == 1][:positive_target]
    negative = [index for index in global_rank if labels[index] == 0][:negative_target]
    indices = positive + negative
    if len(indices) != args.n_ism_sequences:
        raise RuntimeError(
            f"Could select only {len(indices)} of {args.n_ism_sequences} class-balanced ISM examples"
        )
    ism = saturation_mutagenesis(
        model,
        X[indices],
        target=0,
        batch_size=args.ism_batch_size,
        device=device,
        fast=True,
    ).cpu()
    rows = []
    for local_index, example_index in enumerate(indices):
        observed_attr = attribution_scores[example_index].numpy()
        observed_ism = ism[local_index].sum(dim=0).numpy()
        correlation = np.nan
        if np.std(observed_attr) > 0 and np.std(observed_ism) > 0:
            correlation = float(np.corrcoef(observed_attr, observed_ism)[0, 1])
        rows.append(
            {
                "example_idx": example_index,
                "pearson_position_correlation": correlation,
                "mean_absolute_attribution": float(np.abs(observed_attr).mean()),
                "mean_absolute_ism": float(np.abs(observed_ism).mean()),
            }
        )
    return ism, pd.DataFrame(rows), indices


def write_json(path: Path, payload: dict[str, Any]) -> None:
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    validate_args(args)
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)
    torch.cuda.manual_seed_all(args.seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    set_seed(args.seed)
    device = torch.device("cuda")
    started = time.time()

    temporary_parent = args.output_dir.parent
    temporary_parent.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=temporary_parent))
    try:
        model, checkpoint = load_model(args.checkpoint, args.weights, device)
        config = checkpoint["config"]
        sequence_length = int(config["sequence_length"])
        X, manifest = sample_and_encode(
            args.positive_fasta,
            args.negative_fasta,
            args.n_per_class,
            sequence_length,
            args.seed,
            args.all_sequences,
        )
        logits, parity = prediction_parity(model, X, device, args.batch_size)
        probabilities = torch.sigmoid(logits).numpy()
        labels = manifest["label"].to_numpy()
        manifest["logit"] = logits.numpy()
        manifest["probability"] = probabilities
        manifest.to_csv(work / "manifest.tsv", sep="\t", index=False)

        method, audit = run_attribution_audit(model, X, args, device)
        audit["prediction_parity"] = parity
        write_json(work / "audit.json", audit)

        attribution, references, attribution_warnings = attribute_sequences(
            model, X, method, args, device, args.seed
        )
        if method == "deep_lift_shap":
            full_convergence = convergence_metrics(
                model, X, attribution, references, device
            )
            full_convergence["passed"] = (
                full_convergence["relative_error"] <= args.deeplift_tolerance
            )
            audit["full_panel_deeplift"] = full_convergence
            if not full_convergence["passed"]:
                if args.method == "deep_lift_shap":
                    raise RuntimeError(
                        "Full-panel DeepLIFT convergence failed: "
                        f"{full_convergence}"
                    )
                method = "saliency"
                audit["selected_method"] = method
                audit["fallback_reason"] = "full_panel_deeplift_convergence"
                attribution, references, attribution_warnings = attribute_sequences(
                    model, X, method, args, device, args.seed
                )
        write_json(work / "audit.json", audit)
        torch.save(attribution, work / "attributions.pt")
        if references is not None:
            torch.save(references, work / "references.pt")

        attribution_scores = attribution.sum(dim=1)
        screened_candidates = select_candidates(
            X, attribution_scores, manifest, args
        )
        motif_names = list(read_meme(str(args.motif_meme)).keys())
        cluster_metadata = read_cluster_metadata(args.motif_metadata)
        all_seqlets = annotate_windows(
            X,
            screened_candidates,
            args.motif_meme,
            motif_names,
            args.tomtom_pvalue,
            args.n_jobs,
        )
        identity = manifest[["example_idx", "sequence_id", "label", "class_name"]]
        all_seqlets = all_seqlets.merge(
            identity,
            on="example_idx",
            how="left",
        )
        annotation_counts_by_sequence, motif_annotation_summary = (
            summarize_all_annotations(
                all_seqlets,
                manifest,
                motif_names,
                cluster_metadata,
            )
        )

        # The vignette-style motif inventory above uses every called seqlet.
        # Only the expensive supervisor-style ablation is capped and balanced.
        candidates = cap_candidates_stratified(
            all_seqlets, manifest, args.max_candidate_sites
        )
        candidates.insert(0, "candidate_id", np.arange(len(candidates), dtype=int))
        tested_candidates_input = (
            candidates.loc[candidates["motif"].notna()].copy()
            if args.ablate_annotated_only
            else candidates.copy()
        )
        if tested_candidates_input.empty:
            raise RuntimeError(
                "No candidate windows remain for ablation after motif annotation"
            )
        controls = select_matched_controls(
            X,
            attribution_scores,
            tested_candidates_input,
            manifest,
            args,
            excluded_candidates=screened_candidates,
        )
        controls = controls.merge(
            identity,
            on="example_idx",
            how="left",
        )
        tested_candidates_input["site_id"] = (
            "candidate:" + tested_candidates_input["candidate_id"].astype(str)
        )
        tested_candidates_input["site_type"] = "candidate"
        controls["site_id"] = "control:" + controls["control_id"].astype(str)
        controls["site_type"] = "control"
        screened_candidates.to_csv(
            work / "screened_candidate_windows.tsv", sep="\t", index=False
        )
        all_seqlets.to_csv(
            work / "all_seqlets_annotated.tsv", sep="\t", index=False
        )
        annotation_counts_by_sequence.to_csv(
            work / "motif_annotation_counts_by_sequence.tsv", sep="\t", index=False
        )
        motif_annotation_summary.to_csv(
            work / "motif_annotation_summary.tsv", sep="\t", index=False
        )
        candidates.to_csv(work / "candidate_windows.tsv", sep="\t", index=False)
        tested_candidates_input.to_csv(
            work / "tested_candidate_windows.tsv", sep="\t", index=False
        )
        controls.to_csv(work / "matched_control_windows.tsv", sep="\t", index=False)
        logo_manifest = write_attribution_logos(
            attribution,
            all_seqlets,
            candidates,
            manifest,
            args,
            work,
        )
        logo_manifest.to_csv(
            work / "attribution_logo_manifest.tsv", sep="\t", index=False
        )

        all_sites = pd.concat(
            (tested_candidates_input, controls), ignore_index=True, sort=False
        )
        ablation_sites, ablation_replicates = run_ablation(
            model, X, all_sites, args, device
        )
        candidate_ablation = ablation_sites.loc[
            ablation_sites["site_type"] == "candidate"
        ].copy()
        control_ablation = ablation_sites.loc[
            ablation_sites["site_type"] == "control"
        ].copy()
        tested_candidates, control_null = test_candidate_regions(
            candidate_ablation, control_ablation, args
        )
        motif_summary = summarize_clustered_motifs(
            tested_candidates, manifest, cluster_metadata, args
        )
        tested_candidates.to_csv(
            work / "candidate_region_tests.tsv", sep="\t", index=False
        )
        control_ablation.to_csv(
            work / "matched_control_ablation.tsv", sep="\t", index=False
        )
        control_null.to_csv(work / "control_null.tsv", sep="\t", index=False)
        ablation_replicates.to_csv(
            work / "ablation_replicates.tsv", sep="\t", index=False
        )
        motif_summary.to_csv(
            work / "motif_cluster_summary.tsv", sep="\t", index=False
        )

        ism, ism_summary, ism_indices = run_ism_validation(
            model, X, attribution_scores, candidates, manifest, args, device
        )
        torch.save({"example_indices": ism_indices, "attributions": ism}, work / "ism.pt")
        ism_summary = ism_summary.merge(
            manifest[["example_idx", "sequence_id", "label", "class_name"]],
            on="example_idx",
            how="left",
        )
        ism_summary.to_csv(work / "ism_validation.tsv", sep="\t", index=False)

        metadata = {
            "task": args.task,
            "class_semantics": {
                "evn": {"positive": "enhancer", "negative": "genomic background"},
                "svn": {"positive": "silencer", "negative": "genomic background"},
                "evs": {"positive": "enhancer", "negative": "silencer"},
            }[args.task],
            "checkpoint": str(args.checkpoint),
            "weights": str(args.weights),
            "positive_fasta": str(args.positive_fasta),
            "negative_fasta": str(args.negative_fasta),
            "motif_meme": str(args.motif_meme),
            "motif_metadata": str(args.motif_metadata),
            "checksums": {
                "checkpoint": sha256(args.checkpoint),
                "weights": sha256(args.weights),
                "positive_fasta": sha256(args.positive_fasta),
                "negative_fasta": sha256(args.negative_fasta),
                "motif_meme": sha256(args.motif_meme),
                "motif_metadata": sha256(args.motif_metadata),
            },
            "selected_attribution_method": method,
            "parameters": {
                key: value
                for key, value in vars(args).items()
                if key not in {"checkpoint", "weights", "positive_fasta", "negative_fasta", "motif_meme", "motif_metadata", "output_dir"}
            },
            "dataset": {
                "n_sequences": len(X),
                "all_sequences": args.all_sequences,
                "positive_sequences": int((manifest["label"] == 1).sum()),
                "negative_sequences": int((manifest["label"] == 0).sum()),
                "sequence_length": sequence_length,
            },
            "performance_on_sample": {
                "auroc": float(roc_auc_score(labels, probabilities)),
                "average_precision": float(average_precision_score(labels, probabilities)),
            },
            "counts": {
                "screened_candidate_windows": len(screened_candidates),
                "all_seqlets": len(all_seqlets),
                "annotated_all_seqlets": int(all_seqlets["motif"].notna().sum()),
                "motif_clusters_observed": int(
                    all_seqlets.loc[all_seqlets["motif"].notna(), "motif"].nunique()
                ),
                "ablation_candidate_windows": len(candidates),
                "annotated_ablation_candidates": int(candidates["motif"].notna().sum()),
                "tested_candidate_windows": len(tested_candidates_input),
                "important_regions": int(tested_candidates["important_region"].sum()),
                "matched_control_windows": len(controls),
                "enriched_motif_clusters": int(motif_summary["enriched_cluster"].sum()) if len(motif_summary) else 0,
                "ablation_sites": len(ablation_sites),
                "ablation_replicates": len(ablation_replicates),
                "ism_sequences": len(ism_indices),
                "attribution_logos": len(logo_manifest),
            },
            "warnings": {
                "attribution": attribution_warnings,
            },
            "software": {
                "python": sys.version,
                "platform": platform.platform(),
                "torch": torch.__version__,
                "numpy": np.__version__,
                "pandas": pd.__version__,
                "alphagenome_pytorch": getattr(alphagenome_pytorch, "__version__", "unknown"),
                "tangermeme": tangermeme.__version__,
                "cuda": torch.version.cuda,
                "gpu": torch.cuda.get_device_name(0),
            },
            "slurm": {
                key: os.environ.get(key)
                for key in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID")
            },
            "elapsed_seconds": time.time() - started,
        }
        write_json(work / "metadata.json", metadata)
        if args.output_dir.exists():
            shutil.rmtree(args.output_dir)
        work.replace(args.output_dir)
        print(json.dumps(metadata, indent=2, sort_keys=True))
    except Exception:
        shutil.rmtree(work, ignore_errors=True)
        raise


if __name__ == "__main__":
    main()
