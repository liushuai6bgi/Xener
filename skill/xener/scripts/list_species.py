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
from xener import Xener


def main():
    parser = argparse.ArgumentParser(description="List available reference species")
    args = parser.parse_args()

    annor = Xener()
    species = annor.blastdb.keys()

    print(f"Available species ({len(species)}):")
    for s in species:
        print(f"  - {s}")


if __name__ == "__main__":
    main()