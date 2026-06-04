#!/usr/bin/env python3
"""List available reference species from the knowledge graph."""

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