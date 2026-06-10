#!/usr/bin/env python3
"""Refine annotation for one or many clusters into subtypes.

CLI wrapper used by the Xener agent skill. Takes a target cluster ID and two
candidate cell types, then splits the cluster into subtypes using Moran's I
gene filtering and binary/argmax cell assignment.

Skill context: invoked after suggest_refine.py during
references/workflows/refinement.md. The --celltype values must come EXACTLY
from celltype_weight.csv for the target cluster - do not invent names.
Output is written to a new adata.obs column (default: xener_refine).

Two ways to run
---------------
1. Single cluster (legacy, unchanged):

     python scripts/refine_cluster.py --input data.h5ad \\
         --markers output/gene_homolo_weight.csv \\
         --cluster-key leiden --cluster-id 3 \\
         --celltype "lateral root cap,columella root cap cell" \\
         --organ Root --outdir output/

   Writes output/refined_3.csv (+ gexf), as before.

2. Batch / complete-annotation (NEW - shares fixed cost across clusters):

     python scripts/refine_cluster.py --input data.h5ad \\
         --markers output/gene_homolo_weight.csv \\
         --cluster-key leiden --plan output/refine_plan.tsv \\
         --organ Root --outdir output/ \\
         --merge-into output/edf_annotation.csv

   The h5ad (multi-GB), the markers table (tens of MB) and the Xener KG
   connection are loaded ONCE and reused for every cluster in the plan,
   instead of being re-paid by a fresh subprocess per cluster. On a 34k-cell
   Arabidopsis root atlas this fixed cost is ~8 s/cluster; amortising it over
   33 clusters removes ~44% of the batch wall-clock (measured: ~10 min -> ~6 min).

   With --merge-into, the per-cell subtypes are written straight into the
   `xener_refine` column of the consolidated annotation CSV and NO per-cluster
   refined_<id>.csv files are produced (pass --keep-per-cluster to keep them).

Plan file format (--plan)
-------------------------
TSV (one job per line):   <cluster_id>\\t<celltype1,celltype2>
JSON:                     [{"cluster": 3, "celltype": ["a","b"]}, ...]
                    or    {"3": ["a","b"], "13": ["c","d"]}
"""

import argparse
import json
import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from _xener_init import build_xener, add_init_config_arg


# --------------------------------------------------------------------- #
# Plan parsing                                                          #
# --------------------------------------------------------------------- #

def read_plan(path):
    """Parse a batch plan file into a list of (cluster_id:str, celltypes:list[str])."""
    p = str(path)
    if p.lower().endswith(".json"):
        data = json.load(open(p, encoding="utf-8"))
        jobs = []
        if isinstance(data, dict):
            for k, v in data.items():
                cts = v if isinstance(v, list) else str(v).split(",")
                jobs.append((str(k), [str(c).strip() for c in cts if str(c).strip()]))
        else:
            for item in data:
                cid = str(item["cluster"])
                v = item["celltype"]
                cts = v if isinstance(v, list) else str(v).split(",")
                jobs.append((cid, [str(c).strip() for c in cts if str(c).strip()]))
        return jobs

    # TSV: "<cluster_id>\t<celltype1,celltype2>"
    jobs = []
    for raw in open(p, encoding="utf-8"):
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            raise ValueError(
                f"Bad plan line (need '<cluster_id>\\t<celltypes>'): {line!r}"
            )
        cid = parts[0].strip()
        cts = [c.strip() for c in parts[1].split(",") if c.strip()]
        jobs.append((cid, cts))
    return jobs


# --------------------------------------------------------------------- #
# Per-cluster work (shared adata / markers / Xener instance)            #
# --------------------------------------------------------------------- #

def warn_if_truncated(markers_df, cluster_id):
    """Loudly warn if --markers looks like the truncated topk_markers.csv.

    Refinement needs the FULL per-cluster homolog set; a thin set collapses the
    split to a single label (a silent, exit-0 failure).
    """
    if "group" in markers_df.columns and str(cluster_id) in \
            markers_df["group"].astype(str).unique():
        n_homolo = (
            markers_df[markers_df["group"].astype(str) == str(cluster_id)]
            ["homolo"].nunique()
        )
        if n_homolo < 50:
            print(
                f"[WARN] cluster {cluster_id}: only {n_homolo} unique homologs in "
                f"--markers. This looks like topk_markers.csv (truncated). Refinement "
                f"works best on the FULL gene_homolo_weight.csv from Step 3 - a thin "
                f"homolog set frequently yields a single-label, no-split result."
            )


def refine_one(annor, adata, markers_df, cluster_id, celltypes, args):
    """Run a single cluster refinement against already-loaded shared resources."""
    warn_if_truncated(markers_df, cluster_id)
    geneCount, diffgeneCount, annotation, gene2celltype_g = annor.refine_single_cluster(
        adata,
        markers_df,
        cluster_key=args.cluster_key,
        cluster_id=cluster_id,
        candidate_celltype=celltypes,
        key_added=args.key_added,
        organ=args.organ,
        moranI_threshold=args.moran_i,
        split_method=args.split_method,
        markergene_method=args.markergene_method,
    )
    return geneCount, diffgeneCount, annotation, gene2celltype_g


def write_per_cluster(cluster_id, annotation, gene2celltype_g, outdir, key_added):
    """Legacy artifact: refined_<id>.csv + refined_<id>_gene2celltype.gexf."""
    output_csv = os.path.join(outdir, f"refined_{cluster_id}.csv")
    annotation.to_csv(output_csv)
    try:
        import networkx as nx
        output_gexf = os.path.join(outdir, f"refined_{cluster_id}_gene2celltype.gexf")
        nx.write_gexf(gene2celltype_g, output_gexf)
        print(f"Gene->celltype graph saved to {output_gexf}")
    except Exception as e:
        print(f"[WARN] Failed to save gene2celltype gexf for cluster {cluster_id}: {e}")
    print(f"Refinement saved to {output_csv}")


def merge_into_annotation(annotation_csv, collected, key_added):
    """Merge per-cluster subtype labels into one consolidated annotation CSV.

    The `key_added` column (default ``xener_refine``) holds the refinement
    result and ONLY the refinement result: a cell gets a value here iff its
    cluster was actually refined. Cells in unrefined clusters are left EMPTY
    (NaN), not back-filled from `xener`.

    Why empty instead of copying `xener`: `xener_refine` answers a different
    question from `xener`. `xener` is the cluster-level annotation for every
    cell; `xener_refine` is "what did sub-cluster splitting decide for this
    cell". Copying the cluster label into unrefined rows conflates the two and
    silently overstates how much of the dataset was refined - a reader can no
    longer tell which 31k cells were split from which 3k were left alone. An
    empty cell is the honest representation of "not refined"; downstream code
    that wants a fully-populated label should coalesce explicitly
    (e.g. ``adata.obs['xener_refine'].fillna(adata.obs['xener'])``), which keeps
    the choice visible rather than baked in.
    """
    import numpy as np

    ann = pd.read_csv(annotation_csv, index_col=0)
    # Start the column as all-NaN; only refined cells receive a value.
    refine_col = pd.Series(np.nan, index=ann.index, dtype=object)

    n_cells = 0
    n_clusters = 0
    for cid, annotation in collected:
        if key_added not in annotation.columns:
            continue
        common = ann.index.intersection(annotation.index)
        refine_col.loc[common] = annotation.loc[common, key_added].astype(str)
        n_cells += len(common)
        n_clusters += 1
    ann[key_added] = refine_col
    ann.to_csv(annotation_csv)
    n_unrefined = int(ann[key_added].isna().sum())
    print(f"Merged {n_clusters} refinement(s) into '{key_added}' of "
          f"{annotation_csv} ({n_cells} cells labeled, "
          f"{n_unrefined} cells left empty = not refined).")


# --------------------------------------------------------------------- #
# CLI                                                                   #
# --------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[1],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument("--markers", required=True,
                        help="Path to gene_homolo_weight.csv (Step 3 output). "
                             "Do NOT pass topk_markers.csv: refinement needs the FULL "
                             "per-cluster homolog set to query the KG for candidate-celltype "
                             "markers. topk_markers.csv is truncated to top_num genes/cluster, "
                             "which starves the split (argmax dumps all cells into one "
                             "candidate -> no sub-population). The function consumes only the "
                             "group/gene/homolo columns; weights are ignored.")
    parser.add_argument("--cluster-key", default="leiden", help="Cluster column name")

    # Single-cluster mode (legacy)
    parser.add_argument("--cluster-id", default=None,
                        help="Cluster ID to refine (single-cluster mode). "
                             "Mutually exclusive with --plan.")
    parser.add_argument("--celltype", default=None,
                        help="Candidate cell types (comma-separated) for --cluster-id.")

    # Batch mode (new)
    parser.add_argument("--plan", default=None,
                        help="Batch mode: path to a TSV/JSON plan mapping each "
                             "cluster_id to its candidate cell types. All clusters are "
                             "refined in ONE process, sharing the loaded h5ad, markers "
                             "table and KG connection. Mutually exclusive with --cluster-id.")
    parser.add_argument("--merge-into", default=None,
                        help="Path to a consolidated annotation CSV (e.g. "
                             "<dataset>_annotation.csv from run_pipeline.py). Refined "
                             "subtypes are written straight into its `--key-added` column "
                             "and per-cluster refined_<id>.csv files are suppressed.")
    parser.add_argument("--keep-per-cluster", action="store_true",
                        help="Also write per-cluster refined_<id>.csv (+ gexf) even when "
                             "--merge-into is used. Off by default.")

    parser.add_argument("--organ", default=None, help="Organ name")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--key-added", default="xener_refine",
                        help="Column name for refined annotation")
    parser.add_argument("--moran-i", type=float, default=0.5,
                        help="Moran's I threshold for gene filtering [-1, 1]")
    parser.add_argument("--split-method", default="argmax", choices=["bindiv", "argmax"],
                        help="Cluster split method (default argmax: highest refinement success rate)")
    parser.add_argument("--markergene-method", default="all", choices=["diff", "all"],
                        help="Marker gene method (default all: highest refinement success rate)")
    parser.add_argument("--strict", type=int, default=0,
                        help="Strict mode: 0=default, >0=keep max-confidence cell type per gene")
    add_init_config_arg(parser)
    args = parser.parse_args()

    # ---- validate mode selection ----
    if bool(args.plan) == bool(args.cluster_id):
        parser.error("provide EITHER --plan (batch) OR --cluster-id + --celltype "
                     "(single), not both and not neither.")
    if args.cluster_id and not args.celltype:
        parser.error("--celltype is required with --cluster-id.")

    os.makedirs(args.outdir, exist_ok=True)

    # ---- build the job list ----
    if args.plan:
        jobs = read_plan(args.plan)
        if not jobs:
            parser.error(f"plan file {args.plan} contained no jobs.")
        print(f"Batch mode: {len(jobs)} cluster(s) to refine in a single process.")
    else:
        jobs = [(str(args.cluster_id), args.celltype.split(","))]

    # ---- pay the fixed cost ONCE: load h5ad, markers, ensure embedding, connect KG ----
    print(f"Loading h5ad: {args.input}")
    adata = sc.read(args.input)
    print(f"Loading markers: {args.markers}")
    markers_df = pd.read_csv(args.markers)

    if "pca" not in adata.uns:
        print("Running PCA...")
        sc.pp.pca(adata)
    if "connectivities" not in adata.obsp:
        print("Computing neighbors...")
        sc.pp.neighbors(adata)

    annor = build_xener(args.init_config)  # connects to the KG once, reused for every cluster

    # ---- refine every job against the shared resources ----
    collected = []      # list of (cluster_id, annotation_df) for merging
    failures = []
    for i, (cid, cts) in enumerate(jobs, 1):
        tag = f"[{i}/{len(jobs)}] cluster {cid}"
        try:
            print(f"=== {tag}: {','.join(cts)} ===")
            geneCount, diffgeneCount, annotation, g = refine_one(
                annor, adata, markers_df, cid, cts, args
            )
            dist = annotation[args.key_added].value_counts().to_dict()
            print(f"{tag} -> {dist}")
            collected.append((cid, annotation))

            # Per-cluster artifacts: only when explicitly wanted, or in legacy
            # single-cluster mode without --merge-into.
            write_files = args.keep_per_cluster or (args.merge_into is None)
            if write_files:
                write_per_cluster(cid, annotation, g, args.outdir, args.key_added)
        except Exception as e:
            print(f"[ERROR] {tag} failed: {e}")
            failures.append(cid)

    # ---- consolidate ----
    if args.merge_into:
        merge_into_annotation(args.merge_into, collected, args.key_added)

    # ---- summary ----
    n_ok = len(collected)
    n_split = sum(1 for _, a in collected if a[args.key_added].nunique() > 1)
    print(f"\nDone: {n_ok}/{len(jobs)} cluster(s) refined "
          f"({n_split} split into >1 subtype, {n_ok - n_split} uniform).")
    if failures:
        print(f"[WARN] {len(failures)} cluster(s) failed: {failures}")
    print(f"Annotation column: {args.key_added}")


if __name__ == "__main__":
    main()
