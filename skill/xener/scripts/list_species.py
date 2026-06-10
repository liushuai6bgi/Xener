#!/usr/bin/env python3
"""List available reference species from the knowledge graph.

CLI wrapper used by the Xener agent skill. Prints the valid model_species
values that can be used in config.yaml.

Skill context: invoked during references/workflows/config-validation.md to
verify a user-supplied `model_species` list before pipeline execution.
Output must be presented to the user in 6-column tables grouped by Plants /
Animals / etc. — see the formatting rule in skill/xener/SKILL.md.
"""

import argparse
from _xener_init import build_xener, add_init_config_arg


def main():
    parser = argparse.ArgumentParser(description="List available reference species")
    add_init_config_arg(parser)
    args = parser.parse_args()

    annor = build_xener(args.init_config)
    species = annor.blastdb.keys()

    print(f"Available species ({len(species)}):")
    for s in species:
        print(f"  - {s}")


if __name__ == "__main__":
    main()