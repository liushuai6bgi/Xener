#!/usr/bin/env python3
"""List available reference species from the BLAST database.

CLI wrapper used by the Xener agent skill. Prints the valid model_species
values that can be used in config.yaml.

Skill context: invoked during references/workflows/config-validation.md to
verify a user-supplied `model_species` list before pipeline execution.
Output is a simple bullet list of species names. The agent should present
them grouped by taxonomic clade (Plants / Animals / etc.) when displaying
to the user — the raw script output is a flat alphabetical list.
"""

import argparse
from _xener_init import build_xener, add_init_config_arg


def main():
    parser = argparse.ArgumentParser(description="List available reference species")
    add_init_config_arg(parser)
    args = parser.parse_args()

    annor = build_xener(init_config=args.init_config)
    species = annor.blastdb.keys()

    print(f"Available species ({len(species)}):")
    for s in species:
        print(f"  - {s}")


if __name__ == "__main__":
    main()