#!/usr/bin/env python3
"""Step 4: Select top-k genes per cluster."""

import argparse
import os
import pandas as pd

from xener import Xener


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to gene_homolo_weight.zip")
    parser.add_argument("--top_num", type=int, default=30, help="Number of top genes per cluster")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--multihomolo", action="store_true", default=True,
                        help="Keep multiple homologs per gene")
    parser.add_argument("--no-multihomolo", dest="multihomolo", action="store_false",
                        help="Keep only top homolog per gene")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    homolo_weights = pd.read_csv(args.input)
    annor = Xener()

    topk, debug_topk = annor.get_topk_gene(homolo_weights, top_num=args.top_num, multihomolo=args.multihomolo)

    output_path = os.path.join(args.outdir, "topk_markers.zip")
    topk.to_csv(output_path, index=False)
    print(f"Top-k markers saved to {output_path}")
    print(f"Shape: {topk.shape}")


if __name__ == "__main__":
    main()