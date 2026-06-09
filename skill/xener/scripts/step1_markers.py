#!/usr/bin/env python3
"""Step 1: Extract marker genes from clustered single-cell data.

CLI wrapper used by the Xener agent skill. Reads an .h5ad file with cluster
labels in adata.obs, runs Xener's marker gene identification, and writes a
marker_gene.csv checkpoint to --outdir.

Skill context: this is the first of five steps in the Xener pipeline. It is
invoked by `scripts/run_pipeline.py` (the full pipeline) or called directly
when iterating on parameters. See skill/xener/references/workflows/step-by-step.md
for the full workflow and references/parameters.md for CLI flag details.
"""

import argparse
import os
import pandas as pd
import scanpy as sc

from xener import Xener


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument("--cluster-key", default="leiden", help="Cluster column in adata.obs")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--preprocess", action="store_true", help="Run quality control preprocessing")
    parser.add_argument("--batch-key", default=None, help="Batch column for batch correction")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    adata = sc.read(args.input)
    annor = Xener()

    markers = annor.get_markers(
        adata,
        cluster_key=args.cluster_key,
        preprocess=args.preprocess,
        batch_key=args.batch_key
    )

    output_path = os.path.join(args.outdir, "marker_gene.csv")
    markers.to_csv(output_path, index=False)
    print(f"Marker genes saved to {output_path}")
    print(f"Shape: {markers.shape}")


if __name__ == "__main__":
    main()