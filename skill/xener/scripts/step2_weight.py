#!/usr/bin/env python3
"""Step 2: Calculate gene weights from marker genes.

CLI wrapper used by the Xener agent skill. Reads marker_gene.csv from Step 1
and computes per-gene weights using either product (--method prod) or sum
(--method sum) of log fold-changes.

Skill context: invoked by run_pipeline.py or directly during
references/workflows/step-by-step.md when iterating on weighting strategy.
Writes marker_weight.csv as input to Step 3.
"""

import argparse
import os
import pandas as pd

from xener import Xener


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to marker_gene.csv")
    parser.add_argument("--method", default="prod", choices=["prod", "sum"],
                        help="Weight calculation method")
    parser.add_argument("--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    markers = pd.read_csv(args.input)
    annor = Xener()

    weights, debug_gw = annor.get_gene_weight(markers, marker_weight_method=args.method)

    output_path = os.path.join(args.outdir, "marker_weight.csv")
    weights.to_csv(output_path, index=False)
    print(f"Gene weights saved to {output_path}")
    print(f"Shape: {weights.shape}")


if __name__ == "__main__":
    main()