#!/usr/bin/env python3
"""Run full xener pipeline from config file."""

import argparse
import yaml
import os
import sys

from xener import Xener


def main():
    parser = argparse.ArgumentParser(description="Run xener full pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    args = parser.parse_args()

    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    os.makedirs(config["outdir"], exist_ok=True)

    annor = Xener()
    cluster2celltype, cluster2max, debug_params = annor.run_from_yaml(args.config)

    print("Pipeline complete.")
    print("Cluster annotations:", cluster2celltype)
    print("Cluster max-init-weight cell types:", cluster2max)


if __name__ == "__main__":
    main()