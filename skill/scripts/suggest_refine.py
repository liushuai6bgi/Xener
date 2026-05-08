#!/usr/bin/env python3
"""Analyze celltype_weight.csv and output top-5 cell types per cluster for semantic deduplication."""

import argparse
import os
import pandas as pd
import json


def main():
    parser = argparse.ArgumentParser(
        description="Extract top-5 cell types per cluster for refinement analysis"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to celltype_weight.csv (output from step5)",
    )
    parser.add_argument(
        "--topk",
        type=int,
        default=5,
        help="Number of top cell types to return per cluster. Default: 5",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory for suggestions JSON",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)

    suggestions = {}

    for cluster in sorted(df["cluster"].unique()):
        cluster_df = df[df["cluster"] == cluster].sort_values(
            "init_weight", ascending=False
        )

        if len(cluster_df) < 2:
            continue

        top_n = cluster_df.head(args.topk)

        suggestions[cluster] = {
            "top_celltypes": [
                {
                    "celltype": row["celltype"],
                    "init_weight": round(row["init_weight"], 4),
                }
                for _, row in top_n.iterrows()
            ]
        }
    suggestions = {int(cluster): info for cluster, info in suggestions.items()}
    # Print results
    print(f"\n=== Top-{args.topk} Cell Types per Cluster ===\n")
    for cluster, info in sorted(
        suggestions.items(),
        key=lambda x: x[0],
    ):
        print(f"Cluster {cluster}:")
        for item in info["top_celltypes"]:
            print(f"  {item['celltype']}: {item['init_weight']}")
        print()

    # Save to JSON if outdir provided
    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)
        output_path = os.path.join(args.outdir, "top_celltypes.json")
        with open(output_path, "w") as f:
            json.dump(suggestions, f, indent=2)
        print(f"\nTop cell types saved to {output_path}")


if __name__ == "__main__":
    main()
