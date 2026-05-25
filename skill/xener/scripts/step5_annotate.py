#!/usr/bin/env python3
"""Step 5: Annotate cell types from top-k markers."""

import argparse
import os
from pathlib import Path
import pandas as pd
import json

from xener import Xener


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to topk_markers.zip")
    parser.add_argument("--outdir", required=True, help="Output directory for annotations")
    parser.add_argument("--organ", default=None, help="Organ name (e.g., leaf, root)")
    parser.add_argument("--threshold", type=int, default=None,
                        help="Cell type confidence threshold")
    parser.add_argument("--candidate-annotation", nargs="*", default=None,
                        help="Restrict cell types to these candidates")
    parser.add_argument("--mode", default="path", choices=["node", "path"],
                        help="Annotation mode: node (single type) or path (developmental trajectory)")
    parser.add_argument("--decay-factor", type=float, default=0.7,
                        help="Weight decay factor for graph propagation")
    parser.add_argument("--ann-strict", type=int, default=0,
                        help="Strict mode for annotation: <0=loose (binarize), 0=default, 1=max per marker, 2=global max")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    topk = pd.read_csv(args.input)
    annor = Xener()

    cluster2celltype, cluster2max, celltype_weight, debug_ann = \
        annor.cell_annotation(
            topk,
            outdir=outdir,
            organ=args.organ,
            threshold=args.threshold,
            mode=args.mode,
            decay_factor=args.decay_factor,
            ann_strict=args.ann_strict
        )

    # Save results
    with open(os.path.join(args.outdir, "cluster2celltype.json"), "w") as f:
        json.dump({str(k): v for k, v in cluster2celltype.items()}, f, indent=2)

    with open(os.path.join(args.outdir, "cluster2max.json"), "w") as f:
        json.dump({str(k): v for k, v in cluster2max.items()}, f, indent=2)

    celltype_weight.to_csv(os.path.join(args.outdir, "celltype_weight.zip"), index=False)

    print("Annotation complete.")
    print("Cluster annotations:", cluster2celltype)
    print("Cluster max-init-weight:", cluster2max)


if __name__ == "__main__":
    main()