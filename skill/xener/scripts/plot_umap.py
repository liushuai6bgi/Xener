#!/usr/bin/env python3
"""Plot UMAP visualization of Xener annotation and refinement results.

CLI wrapper used by the Xener agent skill. Reads an h5ad file (with cluster
labels in adata.obs, Xener annotation in adata.obs['xener'], optional
refinement in adata.obs['xener_refine'], and UMAP coordinates in
adata.obsm['X_umap']) and writes PNG figures.

Skill context: invoked after Step 5 (or after refinement) to produce
publication-ready UMAP plots. Two modes:

  --mode annotation
      Writes umap_annotation.png with side-by-side panels:
        (left)  cells colored by cluster_key
        (right) cells colored by Xener annotation

  --mode refine --cluster-id N
      Writes umap_refine_cluster_N.png with three panels:
        (left)   all cells colored by cluster_key, cluster N highlighted
        (middle) all cells colored by xener annotation, cluster N highlighted
        (right)  cells in cluster N colored by xener_refine subtype

If --refine-key is provided, the right panel uses that column instead of
xener_refine. If the column does not exist, the right panel is omitted
with a printed warning.

Usage examples:
    # Plot Xener annotation
    python scripts/plot_umap.py --input data.h5ad --mode annotation \\
        --cluster-key leiden --outdir output/

    # Plot refinement for cluster 4
    python scripts/plot_umap.py --input data.h5ad --mode refine \\
        --cluster-key leiden --cluster-id 4 \\
        --refine-key xener_refine --outdir output/

Outputs:
    umap_annotation.png             (--mode annotation)
    umap_refine_cluster_<N>.png     (--mode refine)
"""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # non-interactive backend, safe for headless runs
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns


# --------------------------------------------------------------------- #
# Color helpers                                                         #
# --------------------------------------------------------------------- #

def _palette_for(categories):
    """Pick a perceptually distinct palette scaled to the number of categories."""
    n = len(categories)
    if n <= 10:
        return sns.color_palette("bright", n)
    if n <= 20:
        return sns.color_palette("hls", n)
    return sns.color_palette("husl", n)


def _point_size(n_cells):
    """Pick a point size that avoids heavy overlap."""
    if n_cells < 1000:
        return 50
    if n_cells < 5000:
        return 30
    if n_cells < 10000:
        return 20
    return 10


# --------------------------------------------------------------------- #
# Core scatter                                                          #
# --------------------------------------------------------------------- #

def _scatter(ax, embedding, labels, title, legend=True):
    """Scatter UMAP colored by `labels` (Series aligned with embedding rows)."""
    embedding = np.asarray(embedding)
    if embedding.ndim != 2 or embedding.shape[1] < 2:
        raise ValueError(
            f"UMAP coordinates must be 2D; got shape {embedding.shape}. "
            "Run sc.tl.umap() first or set --embedding-key."
        )

    # Treat labels as strings to avoid categorical→legend weirdness
    labels = labels.astype(str)
    categories = labels.unique().tolist()
    palette = _palette_for(categories)
    color_map = dict(zip(categories, palette))

    s = _point_size(len(labels))
    for cat in categories:
        mask = labels == cat
        ax.scatter(
            embedding[mask, 0],
            embedding[mask, 1],
            c=[color_map[cat]],
            s=s,
            alpha=0.75,
            edgecolors="none",
            label=cat,
        )

    ax.set_xlabel("UMAP_1", fontsize=12)
    ax.set_ylabel("UMAP_2", fontsize=12)
    ax.set_title(title, fontsize=14, pad=12)
    ax.grid(True, alpha=0.25)
    ax.set_aspect("equal", adjustable="datalim")

    if legend:
        # Cap legend to 30 entries to avoid huge legends
        if len(categories) > 30:
            ax.legend([], [], frameon=False)
        else:
            ax.legend(
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                fontsize=8,
                title_fontsize=9,
                frameon=False,
            )


def _scatter_highlight(ax, embedding, labels, highlight_mask, title, legend=True):
    """Like _scatter but only color `highlight_mask` cells; others are gray."""
    embedding = np.asarray(embedding)
    labels = labels.astype(str)

    # Gray background
    ax.scatter(
        embedding[~highlight_mask, 0],
        embedding[~highlight_mask, 1],
        c="lightgray",
        s=_point_size((~highlight_mask).sum()),
        alpha=0.3,
        edgecolors="none",
    )

    # Highlighted cells, colored by label
    sub_labels = labels[highlight_mask]
    categories = sub_labels.unique().tolist()
    palette = _palette_for(categories)
    color_map = dict(zip(categories, palette))
    s = _point_size(highlight_mask.sum())

    for cat in categories:
        idx = np.where(highlight_mask)[0]
        mask_in_sub = sub_labels == cat
        ax.scatter(
            embedding[idx[mask_in_sub], 0],
            embedding[idx[mask_in_sub], 1],
            c=[color_map[cat]],
            s=s,
            alpha=0.85,
            edgecolors="none",
            label=cat,
        )

    ax.set_xlabel("UMAP_1", fontsize=12)
    ax.set_ylabel("UMAP_2", fontsize=12)
    ax.set_title(title, fontsize=14, pad=12)
    ax.grid(True, alpha=0.25)
    ax.set_aspect("equal", adjustable="datalim")

    if legend and len(categories) <= 30:
        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            fontsize=8,
            title_fontsize=9,
            frameon=False,
        )


# --------------------------------------------------------------------- #
# Modes                                                                 #
# --------------------------------------------------------------------- #

def plot_annotation(adata, cluster_key, embedding_key, outdir):
    """Side-by-side: cluster | xener annotation."""
    embedding = adata.obsm[embedding_key]
    cluster_labels = adata.obs[cluster_key].astype(str)
    annot_labels = adata.obs["xener"].astype(str)

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    _scatter(axes[0], embedding, cluster_labels, f"Cluster ({cluster_key})")
    _scatter(axes[1], embedding, annot_labels, "Xener annotation")

    fig.suptitle("Xener Annotation UMAP", fontsize=16, y=1.00)
    fig.tight_layout()
    out_path = Path(outdir) / "umap_annotation.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_refinement(adata, cluster_key, cluster_id, refine_key, embedding_key, outdir):
    """Three-panel: cluster highlight | xener highlight | refine within cluster."""
    embedding = adata.obsm[embedding_key]
    cluster_labels = adata.obs[cluster_key].astype(str)
    annot_labels = adata.obs["xener"].astype(str)

    # Highlight mask: cells whose cluster label matches cluster_id
    cluster_id_str = str(cluster_id)
    highlight = (cluster_labels == cluster_id_str).values
    n_in_cluster = int(highlight.sum())
    if n_in_cluster == 0:
        # Try matching the cluster as a number
        try:
            highlight = (cluster_labels == str(int(cluster_id))).values
            n_in_cluster = int(highlight.sum())
        except ValueError:
            pass
    if n_in_cluster == 0:
        raise ValueError(
            f"No cells found in cluster {cluster_id!r} (column '{cluster_key}'). "
            f"Available: {sorted(cluster_labels.unique())[:10]}..."
        )

    has_refine = refine_key in adata.obs.columns
    if not has_refine:
        print(
            f"[WARN] Column '{refine_key}' not found in adata.obs. "
            "Refinement panel will be skipped. Run refine_cluster.py first."
        )

    n_panels = 3 if has_refine else 2
    fig, axes = plt.subplots(1, n_panels, figsize=(8 * n_panels, 7))

    # Panel 1: cluster view, target cluster highlighted
    _scatter_highlight(
        axes[0], embedding, cluster_labels, highlight,
        f"Cluster {cluster_id_str} (n={n_in_cluster})"
    )

    # Panel 2: xener annotation view, target cluster highlighted
    _scatter_highlight(
        axes[1], embedding, annot_labels, highlight,
        f"Xener annotation (cluster {cluster_id_str} highlighted)"
    )

    if has_refine:
        # Panel 3: only the target cluster, colored by refine
        refine_labels = adata.obs.loc[highlight, refine_key].astype(str)
        _scatter(
            axes[2], embedding[highlight], refine_labels,
            f"Refined subtypes (cluster {cluster_id_str})"
        )

    fig.suptitle(
        f"Refinement of cluster {cluster_id_str}", fontsize=16, y=1.00
    )
    fig.tight_layout()
    out_path = Path(outdir) / f"umap_refine_cluster_{cluster_id_str}.png"
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


# --------------------------------------------------------------------- #
# CLI                                                                   #
# --------------------------------------------------------------------- #

def _load_adata(path, embedding_key):
    """Load the plotting input into an AnnData.

    Accepts either a full ``.h5ad`` or the lightweight annotation ``.csv``
    written by run_pipeline.py (per-cell UMAP coordinates + label columns).
    For the CSV case we construct a minimal AnnData with an empty expression
    matrix and inject the UMAP coordinates into ``obsm`` — no multi-GB matrix
    is ever read or stored.
    """
    if str(path).lower().endswith(".csv"):
        import anndata as ad
        df = pd.read_csv(path, index_col=0)
        adata = ad.AnnData(X=np.zeros((len(df), 1), dtype=np.float32), obs=df.copy())
        if {"UMAP_1", "UMAP_2"}.issubset(df.columns):
            adata.obsm[embedding_key] = df[["UMAP_1", "UMAP_2"]].to_numpy()
        return adata
    return sc.read_h5ad(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument(
        "--mode", required=True, choices=["annotation", "refine"],
        help="annotation: plot Xener result; refine: plot cluster refinement",
    )
    parser.add_argument(
        "--cluster-key", default="leiden",
        help="Column in adata.obs with cluster labels (default: leiden)",
    )
    parser.add_argument(
        "--embedding-key", default="X_umap",
        help="Key in adata.obsm for the 2D embedding (default: X_umap)",
    )
    parser.add_argument(
        "--refine-key", default="xener_refine",
        help="Column with refinement labels (default: xener_refine)",
    )
    parser.add_argument(
        "--cluster-id", default=None,
        help="Cluster to refine (required when --mode refine)",
    )
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"input not found: {args.input}")
    adata = _load_adata(args.input, args.embedding_key)
    if args.embedding_key not in adata.obsm:
        raise KeyError(
            f"Embedding '{args.embedding_key}' not in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )
    if args.cluster_key not in adata.obs:
        raise KeyError(
            f"Cluster key '{args.cluster_key}' not in adata.obs. "
            f"Available: {list(adata.obs.columns)[:10]}"
        )
    if "xener" not in adata.obs:
        raise KeyError(
            "Column 'xener' not found in adata.obs. Run step5_annotate.py first."
        )

    if args.mode == "annotation":
        plot_annotation(adata, args.cluster_key, args.embedding_key, outdir)
    else:
        if args.cluster_id is None:
            parser.error("--cluster-id is required when --mode refine")
        plot_refinement(
            adata, args.cluster_key, args.cluster_id,
            args.refine_key, args.embedding_key, outdir,
        )


if __name__ == "__main__":
    main()
