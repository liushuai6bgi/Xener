#!/usr/bin/env python3
"""List available organs from the knowledge graph."""

import argparse
from xener import Xener


def main():
    parser = argparse.ArgumentParser(description="List available organs")
    parser.add_argument("--species", default=None,
                        help="Filter organs by species (e.g., Brassica_rapa)")
    args = parser.parse_args()

    annor = Xener()
    df = annor.KG.species_organ_cell

    if args.species:
        df = df[df["species"] == args.species]

    organs = sorted(df["organ"].unique())
    print(f"Available organs ({len(organs)}):")
    for o in organs:
        print(f"  - {o}")


if __name__ == "__main__":
    main()