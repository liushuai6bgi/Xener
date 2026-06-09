#!/usr/bin/env python3
"""Refine annotation for a single cluster into subtypes.

CLI wrapper used by the Xener agent skill. Takes a target cluster ID and two
candidate cell types, then splits the cluster into subtypes using Moran's I
gene filtering and binary/argmax cell assignment.

Skill context: invoked after suggest_refine.py during
references/workflows/refinement.md. The --celltype values must come EXACTLY
from celltype_weight.csv for the target cluster — do not invent names.
Output is written to a new adata.obs column (default: xener_refine).
"""

import argparse
import os
import pandas as pd
import scanpy as sc

from xener import Xener

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument("--markers", required=True, help="Path to topk_markers.csv")
    parser.add_argument("--cluster-key", default="leiden", help="Cluster column name")
    parser.add_argument("--cluster-id", required=True, help="Cluster ID to refine")
    parser.add_argument("--celltype", required=True, help="Candidate cell types (comma-separated)")
    parser.add_argument("--organ", default=None, help="Organ name")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--key-added", default="xener_refine",
                        help="Column name for refined annotation")
    parser.add_argument("--moran-i", type=float, default=0.5,
                        help="Moran's I threshold for gene filtering [-1, 1]")
    parser.add_argument("--split-method", default="bindiv", choices=["bindiv", "argmax"],
                        help="Cluster split method")
    parser.add_argument("--markergene-method", default="diff", choices=["diff", "all"],
                        help="Marker gene method")
    parser.add_argument("--strict", type=int, default=0,
                        help="Strict mode: 0=default, >0=keep max-confidence cell type per gene")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.input)
    topk_markers = pd.read_csv(args.markers)

    # Check if PCA and neighbors exist
    if "pca" not in adata.uns:
        print("Running PCA...")
        sc.pp.pca(adata)
    if "connectivities" not in adata.obsp:
        print("Computing neighbors...")
        sc.pp.neighbors(adata)

    annor = Xener()

    geneCount, diffgeneCount, annotation, gene2celltype_g = annor.refine_single_cluster(
        adata,
        topk_markers,
        cluster_key=args.cluster_key,
        cluster_id=args.cluster_id,
        candidate_celltype=args.celltype.split(","),
        key_added=args.key_added,
        organ=args.organ,
        moranI_threshold=args.moran_i,
        split_method=args.split_method,
        markergene_method=args.markergene_method,
    )

    # Save refined annotation
    output_csv = os.path.join(args.outdir, f"refined_{args.cluster_id}.csv")
    annotation.to_csv(output_csv)

    # Save the gene -> homolo -> celltype graph for downstream inspection
    # (e.g., visualizing which marker genes drove the sub-cluster split)
    try:
        import networkx as nx
        output_gexf = os.path.join(args.outdir, f"refined_{args.cluster_id}_gene2celltype.gexf")
        nx.write_gexf(gene2celltype_g, output_gexf)
        print(f"Gene->celltype graph saved to {output_gexf}")
    except Exception as e:
        print(f"[WARN] Failed to save gene2celltype gexf: {e}")

    print(f"Refinement saved to {output_csv}")
    print(f"Annotation column: {args.key_added}")
    print(f"Cell types found: {annotation[args.key_added].unique().tolist()}")
    print(f"Gene counts per cell type: {geneCount}")
    print(f"Diff gene counts per cell type: {diffgeneCount}")
    print(f"Gene->celltype graph nodes: {gene2celltype_g.number_of_nodes()}, edges: {gene2celltype_g.number_of_edges()}")


if __name__ == "__main__":
    main()