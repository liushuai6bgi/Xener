#!/usr/bin/env python3
"""Initialize and validate a custom Xener configuration.

CLI wrapper used by the Xener agent skill. Builds a ``Xener`` object from a
user-supplied init-config YAML and reports what it connected to: the resolved
Knowledge-Graph endpoint, the BLAST database directory, how many reference
species the BLAST DB exposes, and how many organs the KG knows about.

Skill context: run this FIRST whenever a user wants to point Xener at their
own infrastructure (an on-prem Neo4j KG, a locally-built BLAST database, or a
BLASTP result cache) instead of the public cloud defaults. It is the init-time
analogue of ``list_species.py`` / ``list_organs.py`` (which validate the *run*
config): it proves the init-config is usable before the pipeline spends time on
markers and BLAST. See references/workflows/initialization.md.

Init-config schema (all keys optional):

    KG_url:             bolt://localhost:7687     # or http://host:7474; omit -> public cloud KG
    KG_usr:             neo4j                     # KG username
    KG_pwd:             secret                    # KG password
    blastdb_path:       /data/xener/blastdb/prot  # local BLAST protein DB directory
    blastp_result_path: /data/xener/blastp_cache  # optional BLASTP result cache

Usage:
    # Validate a custom configuration (typically before run_pipeline.py)
    python scripts/init_xener.py --init-config my_xener.yaml

    # With no --init-config, this simply checks the default cloud setup works
    python scripts/init_xener.py
"""

import argparse
import sys

# Sibling import: scripts/ is on sys.path[0] when invoked as `python scripts/...`
from _xener_init import build_xener, load_init_config


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument(
        "--init-config", default=None, metavar="PATH",
        help="Path to an init-config YAML (KG_url / KG_usr / KG_pwd, "
             "blastdb_path, blastp_result_path). Omit to validate the default "
             "public cloud KG + bundled BLAST database.",
    )
    args = parser.parse_args()

    if args.init_config:
        # Parse + warn on unknown keys up front, and echo the (password-masked)
        # resolved configuration so the user can confirm what will be used.
        cfg = load_init_config(args.init_config)
        shown = {k: ("***" if k == "KG_pwd" and v else v) for k, v in cfg.items()}
        print(f"Init-config: {args.init_config}")
        print(f"  Resolved keys: {shown if shown else '(empty — using all defaults)'}")
    else:
        print("Init-config: none (default public cloud KG + bundled BLAST database)")

    print("Initializing Xener...")
    try:
        annor = build_xener(args.init_config)
    except Exception as e:
        print(f"[ERROR] Xener initialization failed: {e}", file=sys.stderr)
        print(
            "\nCommon causes:\n"
            "  - KG_url unreachable (wrong host/port, KG not running, firewall)\n"
            "  - bolt:// KG needs KG_usr / KG_pwd and they are missing or wrong\n"
            "  - blastdb_path does not point at a BLAST protein DB directory\n"
            "See references/workflows/initialization.md and "
            "references/troubleshooting.md.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Report what we actually connected to, so the user can sanity-check it.
    species = list(annor.blastdb.keys())
    n_organs = len(annor.KG.available_organ_set)
    print("\nXener initialized successfully.")
    print(f"  BLAST database directory : {annor.blastdb_path}")
    print(f"  Reference species in DB  : {len(species)}")
    print(f"  Organs known to the KG   : {n_organs}")
    if annor.blastp_result_path:
        print(f"  BLASTP result cache      : {annor.blastp_result_path}")

    # Cheap guardrails: an empty BLAST DB or organ set means the endpoints
    # resolved but carry no usable data — better to flag now than mid-pipeline.
    if not species:
        print("\n[WARN] The BLAST database exposes 0 reference species. Check "
              "blastdb_path — it should be a directory of makeblastdb protein "
              "DBs (one per model species).", file=sys.stderr)
    if n_organs == 0:
        print("\n[WARN] The KG returned 0 organs. Check KG_url / credentials — "
              "the endpoint resolved but holds no species_organ_cell data.",
              file=sys.stderr)
    if not species or n_organs == 0:
        sys.exit(1)

    print("\nReady. Pass this same --init-config to run_pipeline.py (or the "
          "step scripts) to run the pipeline against this configuration.")


if __name__ == "__main__":
    main()
