#!/usr/bin/env python3
"""Annotate existing GkmExplain TF-MoDISco patterns with clustered motifs.

This is deliberately a post-discovery stage: it reads previously generated
TF-MoDISco HDF5 files and never invokes TF-MoDISco itself.  Every
pattern--motif comparison is retained so that both per-pattern and global
multiple-testing corrections are auditable.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import h5py
import logomaker
import matplotlib
import numpy as np
import pandas as pd
from memelite import tomtom
from memelite.io import read_meme

matplotlib.use("Agg")
from matplotlib import pyplot as plt  # noqa: E402


DNA = list("ACGT")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=("evn", "svn", "evs"), required=True)
    parser.add_argument("--positive-modisco", type=Path, required=True)
    parser.add_argument("--negative-modisco", type=Path, required=True)
    parser.add_argument("--motif-meme", type=Path, required=True)
    parser.add_argument("--motif-metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return values
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def numeric_pattern_key(name: str) -> tuple[int, str]:
    try:
        return int(name.rsplit("_", 1)[1]), name
    except (IndexError, ValueError):
        return 10**9, name


def load_patterns(path: Path, class_name: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with h5py.File(path, "r") as handle:
        for sign, root_name in (("positive", "pos_patterns"), ("negative", "neg_patterns")):
            if root_name not in handle:
                raise ValueError(f"{path}: missing HDF5 group {root_name}")
            root = handle[root_name]
            for pattern_name in sorted(root.keys(), key=numeric_pattern_key):
                group = root[pattern_name]
                required = ("sequence", "contrib_scores", "hypothetical_contribs", "seqlets")
                missing = [name for name in required if name not in group]
                if missing:
                    raise ValueError(f"{path}:{root_name}/{pattern_name}: missing {missing}")
                sequence = np.asarray(group["sequence"], dtype=float)
                contribution = np.asarray(group["contrib_scores"], dtype=float)
                hypothetical = np.asarray(group["hypothetical_contribs"], dtype=float)
                if sequence.ndim != 2 or sequence.shape[1] != 4:
                    raise ValueError(
                        f"{path}:{root_name}/{pattern_name}: invalid sequence shape {sequence.shape}"
                    )
                if contribution.shape != sequence.shape or hypothetical.shape != sequence.shape:
                    raise ValueError(
                        f"{path}:{root_name}/{pattern_name}: inconsistent pattern matrices"
                    )
                seqlet_group = group["seqlets"]
                if "n_seqlets" in seqlet_group:
                    n_seqlets = int(np.asarray(seqlet_group["n_seqlets"]).reshape(-1)[0])
                elif "example_idx" in seqlet_group:
                    n_seqlets = int(len(seqlet_group["example_idx"]))
                else:
                    raise ValueError(
                        f"{path}:{root_name}/{pattern_name}: cannot determine seqlet count"
                    )
                records.append(
                    {
                        "pattern_id": f"{class_name}:{sign}:{pattern_name}",
                        "class_name": class_name,
                        "sign": sign,
                        "pattern_name": pattern_name,
                        "pattern_length": int(sequence.shape[0]),
                        "n_seqlets": n_seqlets,
                        "sequence": sequence,
                        "contribution": contribution,
                        "hypothetical": hypothetical,
                    }
                )
    return records


def read_cluster_metadata(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    required = {
        "Representative_Motif",
        "Representative_Motif_Gene_Symbol",
        "Cluster_Size",
        "Primary_Family",
        "Primary_Subfamily",
        "Clustered_Motifs_Gene_Symbols",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{path}: missing metadata columns {missing}")
    if table["Representative_Motif"].duplicated().any():
        raise ValueError(f"{path}: duplicate representative motifs")
    return table


def annotate_patterns(
    patterns: list[dict[str, object]],
    motifs: dict[str, np.ndarray],
    n_jobs: int,
    fdr: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    queries = [np.asarray(pattern["sequence"], dtype=float).T for pattern in patterns]
    motif_names = list(motifs)
    result = tomtom(
        queries,
        list(motifs.values()),
        n_nearest=None,
        n_jobs=n_jobs,
    )
    pvalues, scores, offsets, overlaps, strands = result
    if pvalues.shape != (len(patterns), len(motif_names)):
        raise RuntimeError(f"Unexpected TOMTOM result shape: {pvalues.shape}")

    global_qvalues = bh_adjust(pvalues.reshape(-1)).reshape(pvalues.shape)
    query_qvalues = np.vstack([bh_adjust(row) for row in pvalues])
    rows: list[dict[str, object]] = []
    for query_idx, pattern in enumerate(patterns):
        for motif_idx, motif in enumerate(motif_names):
            rows.append(
                {
                    "pattern_id": pattern["pattern_id"],
                    "class_name": pattern["class_name"],
                    "sign": pattern["sign"],
                    "pattern_name": pattern["pattern_name"],
                    "pattern_length": pattern["pattern_length"],
                    "n_seqlets": pattern["n_seqlets"],
                    "motif_index": motif_idx,
                    "motif": motif,
                    "tomtom_pvalue": float(pvalues[query_idx, motif_idx]),
                    "tomtom_query_qvalue": float(query_qvalues[query_idx, motif_idx]),
                    "tomtom_global_qvalue": float(global_qvalues[query_idx, motif_idx]),
                    "tomtom_score": float(scores[query_idx, motif_idx]),
                    "tomtom_offset": int(offsets[query_idx, motif_idx]),
                    "tomtom_overlap": int(overlaps[query_idx, motif_idx]),
                    "tomtom_strand": str(strands[query_idx, motif_idx]),
                }
            )
    matches = pd.DataFrame(rows)
    matches["significant"] = matches["tomtom_global_qvalue"] <= fdr
    matches = matches.sort_values(
        ["pattern_id", "tomtom_pvalue", "motif"], ignore_index=True
    )
    annotations = matches.groupby("pattern_id", sort=False).first().reset_index()
    annotations = annotations.rename(
        columns={
            "motif": "best_motif",
            "motif_index": "best_motif_index",
            "tomtom_pvalue": "best_tomtom_pvalue",
            "tomtom_query_qvalue": "best_tomtom_query_qvalue",
            "tomtom_global_qvalue": "best_tomtom_global_qvalue",
            "tomtom_score": "best_tomtom_score",
            "tomtom_offset": "best_tomtom_offset",
            "tomtom_overlap": "best_tomtom_overlap",
            "tomtom_strand": "best_tomtom_strand",
            "significant": "best_match_significant",
        }
    )
    significant_counts = (
        matches.loc[matches["significant"]]
        .groupby("pattern_id")
        .size()
        .rename("n_significant_motif_matches")
    )
    annotations = annotations.merge(significant_counts, on="pattern_id", how="left")
    annotations["n_significant_motif_matches"] = (
        annotations["n_significant_motif_matches"].fillna(0).astype(int)
    )
    return matches, annotations


def add_cluster_metadata(
    annotations: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    renamed = metadata.rename(
        columns={
            "Representative_Motif": "best_motif",
            "Representative_Motif_Gene_Symbol": "representative_gene",
            "Cluster_Size": "cluster_size",
            "Primary_Family": "primary_family",
            "Primary_Subfamily": "primary_subfamily",
            "Clustered_Motifs_Gene_Symbols": "clustered_gene_symbols",
        }
    )
    columns = [
        "best_motif",
        "representative_gene",
        "cluster_size",
        "primary_family",
        "primary_subfamily",
        "clustered_gene_symbols",
    ]
    return annotations.merge(renamed[columns], on="best_motif", how="left")


def motif_summary(annotations: pd.DataFrame) -> pd.DataFrame:
    accepted = annotations.loc[annotations["best_match_significant"]].copy()
    if accepted.empty:
        return pd.DataFrame(
            columns=[
                "best_motif", "representative_gene", "cluster_size",
                "primary_family", "primary_subfamily", "class_name", "sign",
                "n_patterns", "total_seqlets", "mean_seqlets_per_pattern",
                "minimum_global_qvalue",
            ]
        )
    return (
        accepted.groupby(
            [
                "best_motif", "representative_gene", "cluster_size",
                "primary_family", "primary_subfamily", "class_name", "sign",
            ],
            dropna=False,
            observed=True,
        )
        .agg(
            n_patterns=("pattern_id", "size"),
            total_seqlets=("n_seqlets", "sum"),
            mean_seqlets_per_pattern=("n_seqlets", "mean"),
            minimum_global_qvalue=("best_tomtom_global_qvalue", "min"),
        )
        .reset_index()
        .sort_values(
            ["class_name", "sign", "total_seqlets", "n_patterns"],
            ascending=[True, True, False, False],
            ignore_index=True,
        )
    )


def information_logo_matrix(sequence: np.ndarray) -> pd.DataFrame:
    denominator = sequence.sum(axis=1, keepdims=True)
    probability = np.divide(
        sequence,
        denominator,
        out=np.full_like(sequence, 0.25, dtype=float),
        where=denominator > 0,
    )
    frame = pd.DataFrame(probability, columns=DNA)
    return logomaker.transform_matrix(
        frame,
        from_type="probability",
        to_type="information",
        background=np.full(4, 0.25),
    )


def write_logos(
    patterns: list[dict[str, object]],
    annotations: pd.DataFrame,
    output_dir: Path,
) -> pd.DataFrame:
    logo_dir = output_dir / "logos"
    logo_dir.mkdir(parents=True, exist_ok=True)
    annotation_lookup = annotations.set_index("pattern_id")
    rows = []
    for pattern in patterns:
        pattern_id = str(pattern["pattern_id"])
        annotation = annotation_lookup.loc[pattern_id]
        sequence = np.asarray(pattern["sequence"], dtype=float)
        contribution = np.asarray(pattern["contribution"], dtype=float)
        fig, axes = plt.subplots(2, 1, figsize=(15, 7.5), sharex=True)
        logomaker.Logo(information_logo_matrix(sequence), ax=axes[0])
        axes[0].set_ylabel("Information (bits)", fontsize=14, labelpad=10)
        axes[0].set_title(
            f"{pattern_id} | n={int(pattern['n_seqlets']):,} seqlets | "
            f"{annotation['best_motif']} | global q={annotation['best_tomtom_global_qvalue']:.2g}",
            fontsize=14,
            pad=12,
        )
        logomaker.Logo(pd.DataFrame(contribution, columns=DNA), ax=axes[1])
        axes[1].axhline(0, color="black", linewidth=0.8)
        axes[1].set_ylabel("Contribution", fontsize=14, labelpad=10)
        axes[1].set_xlabel("Pattern position", fontsize=14, labelpad=10)
        tick_step = 5 if sequence.shape[0] <= 30 else 10
        axes[1].set_xticks(np.arange(0, sequence.shape[0], tick_step))
        for axis in axes:
            axis.tick_params(axis="both", labelsize=12, width=0.8)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
        fig.tight_layout(pad=2.0, h_pad=2.0)
        stem = pattern_id.replace(":", "_")
        png = logo_dir / f"{stem}.png"
        svg = logo_dir / f"{stem}.svg"
        fig.savefig(png, dpi=200, bbox_inches="tight")
        fig.savefig(svg, bbox_inches="tight")
        plt.close(fig)
        rows.append(
            {
                "pattern_id": pattern_id,
                "png": str(png.relative_to(output_dir)),
                "svg": str(svg.relative_to(output_dir)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    if args.n_jobs <= 0:
        raise ValueError("--n-jobs must be positive")
    if not 0 < args.fdr <= 1:
        raise ValueError("--fdr must be in (0, 1]")
    inputs = (
        args.positive_modisco,
        args.negative_modisco,
        args.motif_meme,
        args.motif_metadata,
    )
    for path in inputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")
    output_files = (
        "motif_matches.tsv",
        "pattern_annotations.tsv",
        "motif_summary.tsv",
        "logo_manifest.tsv",
        "metadata.json",
    )
    existing = [args.output_dir / name for name in output_files if (args.output_dir / name).exists()]
    if existing and not args.overwrite:
        raise FileExistsError(f"Output already exists: {existing[0]}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=".gkmexplain_annotation.", dir=args.output_dir))
    try:
        patterns = load_patterns(args.positive_modisco, "positive")
        patterns.extend(load_patterns(args.negative_modisco, "negative"))
        if not patterns:
            raise RuntimeError("No TF-MoDISco patterns found")
        motifs = read_meme(str(args.motif_meme))
        cluster_metadata = read_cluster_metadata(args.motif_metadata)
        missing_metadata = sorted(set(motifs) - set(cluster_metadata["Representative_Motif"]))
        if missing_metadata:
            raise ValueError(
                f"Metadata missing {len(missing_metadata)} motif representatives; "
                f"first: {missing_metadata[0]}"
            )

        matches, annotations = annotate_patterns(
            patterns, motifs, args.n_jobs, args.fdr
        )
        annotations = add_cluster_metadata(annotations, cluster_metadata)
        summary = motif_summary(annotations)
        logos = write_logos(patterns, annotations, work)

        matches.to_csv(work / "motif_matches.tsv", sep="\t", index=False)
        annotations.to_csv(work / "pattern_annotations.tsv", sep="\t", index=False)
        summary.to_csv(work / "motif_summary.tsv", sep="\t", index=False)
        logos.to_csv(work / "logo_manifest.tsv", sep="\t", index=False)
        metadata = {
            "task": args.task,
            "discovery_rerun": False,
            "annotation_only": True,
            "positive_modisco": str(args.positive_modisco),
            "negative_modisco": str(args.negative_modisco),
            "motif_meme": str(args.motif_meme),
            "motif_metadata": str(args.motif_metadata),
            "fdr": args.fdr,
            "multiple_testing": {
                "query_qvalue": "Benjamini-Hochberg within each pattern across all motifs",
                "global_qvalue": "Benjamini-Hochberg across all pattern-motif comparisons",
                "significance_column": "tomtom_global_qvalue",
            },
            "counts": {
                "patterns": len(patterns),
                "motifs": len(motifs),
                "comparisons": len(matches),
                "significant_patterns": int(annotations["best_match_significant"].sum()),
            },
            "checksums": {
                str(path): sha256(path) for path in inputs
            },
        }
        with (work / "metadata.json").open("w") as handle:
            json.dump(metadata, handle, indent=2, sort_keys=True)
            handle.write("\n")

        for name in output_files:
            destination = args.output_dir / name
            if destination.exists():
                destination.unlink()
            (work / name).replace(destination)
        logo_destination = args.output_dir / "logos"
        if logo_destination.exists():
            shutil.rmtree(logo_destination)
        (work / "logos").replace(logo_destination)
        print(json.dumps(metadata, indent=2, sort_keys=True))
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
