#!/usr/bin/env python3
"""Shared Xener-construction helper for the skill's scripts.

NOT a CLI entry point — this module is imported by the other scripts; do not
run it directly. (To test a custom configuration, use ``init_xener.py``.)

Why this exists
---------------
xener has TWO different YAML configs, and they answer two different questions:

  * **Run config** (``config.yaml`` -> ``run_pipeline.py``) says *WHAT* to
    annotate: ``model_species``, ``organ``, ``non_model_h5ad``, ``top_num`` ...
    This is the file every user already edits.

  * **Init config** (this module) says *WHERE* Xener gets its data: which
    Knowledge Graph to connect to and which BLAST database to use. It is what
    the ``Xener(...)`` constructor takes. Until now every script built a bare
    ``Xener()``, so it always used the public cloud KG (https://xenor.dcs.cloud)
    and the bundled BLAST database. A user running an on-prem Neo4j KG, or
    pointing at a locally-built BLAST DB, had no way to say so.

This helper is the single place a ``Xener`` object is built, so EVERY script
in the skill can honor the same optional init config in the same way.

Init-config schema (a small YAML; **all keys optional** — omit a key to keep
its default):

    KG_url:             bolt://localhost:7687     # or http://host:7474, or omit for the public cloud
    KG_usr:             neo4j                     # KG username (bolt backends)
    KG_pwd:             secret                    # KG password
    blastdb_path:       /data/xener/blastdb/prot  # local BLAST protein DB directory
    blastp_result_path: /data/xener/blastp_cache  # optional BLASTP result cache

Backward compatible: with no init config, ``build_xener()`` is exactly
``Xener()`` — the public cloud KG plus the bundled BLAST database, as before.

See references/workflows/initialization.md for the full workflow.
"""

# Init-config keys that may be inlined directly into a run config.yaml
# (run_pipeline.py extracts these when no separate init-config file is given).
INIT_CONFIG_KEYS = ("KG_url", "KG_usr", "KG_pwd", "blastdb_path", "blastp_result_path")


def add_init_config_arg(parser):
    """Register the shared ``--init-config`` flag on an argparse parser.

    Call this in every script that builds a Xener so the flag — and its help
    text — stays identical everywhere.
    """
    parser.add_argument(
        "--init-config", default=None, metavar="PATH",
        help="Path to an init-config YAML describing WHERE Xener gets its data "
             "(KG_url / KG_usr / KG_pwd, blastdb_path, blastp_result_path). "
             "Omit to use the public cloud KG + bundled BLAST database. "
             "See references/workflows/initialization.md.",
    )


def load_init_config(path):
    """Read and lightly validate an init-config YAML; return a dict.

    Decoded as UTF-8 explicitly (see mandatory-rules.md §2b) so a non-ASCII
    comment never trips the platform-default codec on a zh-CN Windows console.
    """
    import yaml

    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}
    if not isinstance(cfg, dict):
        raise ValueError(
            f"init-config {path} must be a YAML mapping (key: value pairs); "
            f"got {type(cfg).__name__}."
        )
    unknown = set(cfg) - set(INIT_CONFIG_KEYS)
    if unknown:
        # Not fatal — a user may keep notes in the file — but surface it, since
        # a typo'd key (e.g. `KG_URL`) would otherwise be silently ignored.
        print(f"[WARN] init-config {path}: ignoring unrecognized key(s) "
              f"{sorted(unknown)}. Known keys: {list(INIT_CONFIG_KEYS)}.")
    return cfg


def build_xener(init_config=None):
    """Construct a ``Xener`` instance, optionally from a user init-config.

    Parameters
    ----------
    init_config : str | dict | None
        Path to an init-config YAML, an already-parsed dict of the same keys,
        or ``None`` for defaults (public cloud KG + bundled BLAST database).

    Notes
    -----
    Tolerant of partial configs, unlike the package's ``Xener.init_from_yaml``
    (which requires both ``KG_url`` and ``blastdb_path``). Any subset of keys
    works: a config that sets only ``KG_url`` keeps the bundled BLAST DB; one
    that sets only ``blastdb_path`` keeps the cloud KG. A missing or empty key
    falls back to the same default the bare ``Xener()`` constructor would use.
    """
    from xener import Xener

    if init_config is None:
        return Xener()

    cfg = init_config if isinstance(init_config, dict) else load_init_config(init_config)

    # Only build kg_kwargs when the user actually specified a KG endpoint;
    # otherwise leave it None so Xener falls back to its default cloud url.
    kg_kwargs = None
    if cfg.get("KG_url"):
        kg_kwargs = {
            "url": cfg["KG_url"],
            "usr": cfg.get("KG_usr"),
            "pwd": cfg.get("KG_pwd"),
        }

    return Xener(
        kg_kwargs=kg_kwargs,
        blastdb_path=cfg.get("blastdb_path"),
        blastp_result_path=cfg.get("blastp_result_path"),
    )
