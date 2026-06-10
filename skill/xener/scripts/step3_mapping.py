#!/usr/bin/env python3
"""Step 3: BLAST homology mapping across species.

CLI wrapper used by the Xener agent skill. Maps each marker gene to the
chosen model species via BLASTP, filtering by --pident, --evalue, and
--bitscore thresholds.

Skill context: invoked by run_pipeline.py or directly during
references/workflows/step-by-step.md. Writes blastp_{species}.csv and
gene_homolo_weight.csv as input to Step 4. Pass --species multiple times
(once per model species) to combine homology evidence across references.
"""

import argparse
import os
from pathlib import Path
import pandas as pd

from _xener_init import build_xener, add_init_config_arg


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
    parser.add_argument("--mapping-strict", type=int, default=0,
                        help="Strict mode for BLAST mapping: <0=loose (weight=1), 0=default, 1=top per group/gene")
    add_init_config_arg(parser)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    markers = pd.read_csv(args.input)
    annor = build_xener(args.init_config)

    homolo_weights, debug_map = annor.mapping(
        markers,
        non_model_fasta=args.fasta,
        model_species=args.species,
        outdir=outdir,
        homolo_weight_key=args.weight_key,
        pident=args.pident,
        evalue=args.evalue,
        bitscore=args.bitscore,
        num_threads=args.num_threads,
        mapping_strict=args.mapping_strict
    )

    output_path = os.path.join(args.outdir, "gene_homolo_weight.csv")
    homolo_weights.to_csv(output_path, index=False)
    print(f"Homology mapping saved to {output_path}")
    print(f"Shape: {homolo_weights.shape}")


if __name__ == "__main__":
    main()