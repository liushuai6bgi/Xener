#!/usr/bin/env python3
"""Step 4: Select top-k genes per cluster.

CLI wrapper used by the Xener agent skill. Reads gene_homolo_weight.csv from
Step 3 and selects the top K genes per cluster (default K=30). Set
--multihomolo to keep multiple homologs per non-model gene.

Skill context: invoked by run_pipeline.py or directly during
references/workflows/step-by-step.md when tuning --k. Writes
topk_markers.csv as input to Step 5.
"""

import argparse
import os
import pandas as pd

from _xener_init import build_xener, add_init_config_arg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to gene_homolo_weight.csv")
    parser.add_argument("--top_num", type=int, default=30, help="Number of top genes per cluster")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--multihomolo", action="store_true", default=True,
                        help="Keep multiple homologs per gene")
    parser.add_argument("--no-multihomolo", dest="multihomolo", action="store_false",
                        help="Keep only top homolog per gene")
    add_init_config_arg(parser)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    homolo_weights = pd.read_csv(args.input)
    annor = build_xener(args.init_config)

    topk, debug_topk = annor.get_topk_gene(homolo_weights, top_num=args.top_num, multihomolo=args.multihomolo)

    output_path = os.path.join(args.outdir, "topk_markers.csv")
    topk.to_csv(output_path, index=False)
    print(f"Top-k markers saved to {output_path}")
    print(f"Shape: {topk.shape}")


if __name__ == "__main__":
    main()