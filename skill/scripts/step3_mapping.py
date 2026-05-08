#!/usr/bin/env python3
"""Step 3: BLAST homology mapping across species."""

import argparse
import os
from pathlib import Path
import pandas as pd

from xener import Xener


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to marker_weight.csv")
    parser.add_argument("--fasta", required=True, help="Path to non-model species FASTA file")
    parser.add_argument("--species", required=True, nargs="+",
                        help="Model species names (e.g., Brassica_rapa)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--weight-key", default="pident",
                        help="Column to use as homology weight (pident, evalue, bitscore)")
    parser.add_argument("--pident", type=float, default=60,
                        help="Minimum percent identity filter")
    parser.add_argument("--evalue", type=float, default=0.05,
                        help="Maximum evalue filter")
    parser.add_argument("--bitscore", type=float, default=200,
                        help="Minimum bitscore filter")
    parser.add_argument("--num-threads", type=int, default=None,
                        help="Number of threads for BLAST")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    markers = pd.read_csv(args.input)
    annor = Xener()

    homolo_weights = annor.mapping(
        markers,
        non_model_fasta=args.fasta,
        model_species=args.species,
        outdir=outdir,
        as_homolo_weight_key=args.weight_key,
        pident=args.pident,
        evalue=args.evalue,
        bitscore=args.bitscore,
        num_threads=args.num_threads
    )

    output_path = os.path.join(args.outdir, "gene_homolo_weight.csv")
    homolo_weights.to_csv(output_path, index=False)
    print(f"Homology mapping saved to {output_path}")
    print(f"Shape: {homolo_weights.shape}")


if __name__ == "__main__":
    main()