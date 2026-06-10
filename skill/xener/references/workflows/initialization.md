# Initialization Workflow (custom KG / BLAST database)

By default every script connects Xener to the **public cloud Knowledge Graph**
(`https://xenor.dcs.cloud`) and the **bundled BLAST database** (downloaded to
`~/.xener/data/` on first use). Most users never need anything else — skip this
workflow entirely and the pipeline just works.

Read this only when a user wants Xener to use **their own infrastructure**:
an on-prem Neo4j Knowledge Graph, a locally-built BLAST protein database, or a
BLASTP result cache.

## Two configs, two questions — do not conflate them

Xener is driven by two separate YAML files that answer different questions:

| File | Question | Keys | Passed to |
|------|----------|------|-----------|
| **run config** (`config.yaml`) | *WHAT* do I annotate? | `model_species`, `organ`, `non_model_h5ad`, `non_model_fasta`, `cluster_key`, `top_num`, ... | `--config` |
| **init config** (`xener-init.yaml`) | *WHERE* does Xener get its data? | `KG_url`, `KG_usr`, `KG_pwd`, `blastdb_path`, `blastp_result_path` | `--init-config` |

The run config is documented in `config-schema.md`. This file documents the
**init config**.

## Init-config schema (all keys optional)

```yaml
# Knowledge Graph endpoint. Omit all three to use the public cloud KG.
KG_url: bolt://localhost:7687     # http(s):// -> HTTP backend; bolt:// -> Neo4j Bolt backend
KG_usr: neo4j                     # KG username (Bolt backend only)
KG_pwd: secret                    # KG password (Bolt backend only)

# Local BLAST protein DB directory (makeblastdb output, one DB per model
# species). Omit to use the bundled database.
blastdb_path: /data/xener/blastdb/prot

# Optional BLASTP result cache directory (reuse BLAST results across runs).
blastp_result_path: /data/xener/blastp_cache
```

Every key is optional and independent. A config that sets **only** `KG_url`
keeps the bundled BLAST DB; one that sets **only** `blastdb_path` keeps the
cloud KG. A missing key falls back to the same default a bare `Xener()` would
use. An example file ships at `examples/xener-init.example.yaml`.

> Keep comments ASCII and save as UTF-8 — see `mandatory-rules.md` §2b. The
> init config is read with an explicit UTF-8 decode, like the run config.

## Step 1 — Validate the init config (do this first)

`init_xener.py` is the init-time analogue of `list_species.py`: it builds the
Xener object from your config and reports what it connected to, **before** the
pipeline spends time on markers and BLAST.

```bash
python scripts/init_xener.py --init-config xener-init.yaml
```

A healthy result prints the resolved BLAST directory, the number of reference
species in the DB, and the number of organs the KG knows about, then exits 0:

```
Xener initialized successfully.
  BLAST database directory : /data/xener/blastdb/prot
  Reference species in DB  : 27
  Organs known to the KG   : 8
```

If it exits non-zero, the message names the likely cause (KG unreachable,
missing Bolt credentials, or a `blastdb_path` that is not a BLAST DB). Fix the
config and re-run this step until it passes — see `troubleshooting.md`.

> The reference-species count here is also the menu for `model_species` in your
> run config, and the organ count is the menu for `organ`. Confirm those values
> with `list_species.py --init-config xener-init.yaml` and
> `list_organs.py --init-config xener-init.yaml` (the **same** `--init-config`),
> so validation queries the same KG the pipeline will use.

## Step 2 — Run the pipeline against the custom config

Pass the **same** `--init-config` to the pipeline. Two equivalent ways:

**A. Separate init file (recommended — keeps "what" and "where" apart):**

```bash
python scripts/run_pipeline.py --config config.yaml --init-config xener-init.yaml
```

**B. Inlined into the run config** (convenient for a one-off; just add the
init keys to `config.yaml`):

```yaml
# config.yaml — run keys PLUS init keys in one file
cluster_key: leiden
model_species: [Arabidopsis_thaliana]
organ: Root
non_model_h5ad: edf.h5ad
non_model_fasta: abc.fasta
outdir: output/edf
# --- init keys (WHERE), read by run_pipeline.py ---
KG_url: bolt://localhost:7687
KG_usr: neo4j
KG_pwd: secret
blastdb_path: /data/xener/blastdb/prot
```

```bash
python scripts/run_pipeline.py --config config.yaml
```

A separate `--init-config` file takes precedence over inlined keys when both
are present.

## Using a custom config with the step / refinement scripts

Every script that builds a Xener accepts the same `--init-config` flag:
`list_species.py`, `list_organs.py`, `step1_markers.py` … `step5_annotate.py`,
and `refine_cluster.py`. Pass it consistently across a multi-step run so every
stage talks to the same KG and BLAST DB:

```bash
python scripts/step3_mapping.py --input output/marker_weight.csv \
    --fasta abc.fasta --species Arabidopsis_thaliana \
    --outdir output/ --init-config xener-init.yaml

python scripts/refine_cluster.py --input edf.h5ad \
    --markers output/gene_homolo_weight.csv --cluster-key leiden \
    --plan output/refine_plan.tsv --organ Root --outdir output/ \
    --merge-into output/edf_annotation.csv --init-config xener-init.yaml
```

## Mental model

`scripts/_xener_init.py` is the single place a `Xener` object is built. Every
script imports `build_xener()` from it, so the init config is honored
identically everywhere. With no init config, `build_xener()` is exactly
`Xener()` — fully backward compatible.
