#!/usr/bin/env python3
"""List available organs from the knowledge graph.

CLI wrapper used by the Xener agent skill. Prints the valid `organ` values
that can be used in config.yaml (optionally filtered by --species).

Skill context: invoked during references/workflows/config-validation.md to
verify a user-supplied `organ` value before pipeline execution. Output must
be presented to the user for explicit confirmation.
"""

import argparse
from _xener_init import build_xener, add_init_config_arg


def main():
    parser = argparse.ArgumentParser(description="List available organs")
    parser.add_argument("--species", default=None,
                        help="Filter organs by species (e.g., Brassica_rapa)")
    add_init_config_arg(parser)
    args = parser.parse_args()

    annor = build_xener(args.init_config)
    df = annor.KG.species_organ_cell

    if args.species:
        df = df[df["species"] == args.species]

    organs = sorted(df["organ"].unique())
    print(f"Available organs ({len(organs)}):")
    for o in organs:
        print(f"  - {o}")


if __name__ == "__main__":
    main()