# h5ad Inspection Workflow

Before running the pipeline, the agent **must inspect the h5ad file**
to extract metadata and make informed parameter decisions. This file
documents the inspection procedure.

## Why inspect first?

- The agent needs `cluster_key` — usually stored in `adata.obs`
- The agent may need species hints (var name prefix, uns metadata)
- The agent may need organ hints (uns metadata, file path)
- The agent needs to estimate `top_num` based on cluster sizes

## Inspect with `scripts/inspect_h5ad.py` (load the h5ad ONCE)

Run the dedicated inspection script. It loads the (often multi-GB) h5ad
**exactly once** and prints every field you need to compose a config:

```bash
python scripts/inspect_h5ad.py <path_to_h5ad>
python scripts/inspect_h5ad.py <path_to_h5ad> --json   # machine-readable
```

It reports: shape; obs columns + cardinality; **recommended `cluster_key`**
(leiden → louvain → *cluster*, skipping existing `celltype`/`annotation`
columns); `uns` keys; `var_names` sample **+ detected species**; `obsm` keys
(is `X_umap` present?); whether `raw` is set; an `.X` sparsity/range/log-norm
sanity check; **cluster-size summary** (n / min / median / mean / max) for the
recommended key; a **recommended `top_num`**; and **organ/tissue hints**.

> **Do NOT follow this with ad-hoc `python -c "import scanpy; sc.read(...)"`
> reads to grab "one more field".** Every field needed to write the config is
> already in the output above. Each extra read of a ~1 GB h5ad is ~30-45 s and
> is the single largest I/O cost of the workflow. This is the I/O discipline
> mandated by `mandatory-rules.md` §10. (`inspect_h5ad.py` imports `scanpy`
> directly; this is the sanctioned inspection path, alongside the rule §1
> exception. It never `import xener` and needs no `--init-config`.)

The decision tables below explain how to *interpret* the script's output.

## Decision rules

### cluster_key

Pick the first column with cardinality 5-50 (typical clustering range):

| Column name contains | Priority |
|---------------------|----------|
| `leiden`, `louvain` | High |
| `celltype`, `annotation` (existing) | Skip (don't re-annotate) |
| `cluster`, `cluster_id` | Medium |
| Numeric only (e.g., `X_kmeans`) | Low |

**Fallback**: if no clustering column exists, run leiden first
(agent should add a step before the pipeline).

### Species hint detection

| Signal | Detected species |
|--------|------------------|
| var_names start with `AT` or `AT[0-9]G` | Arabidopsis_thaliana |
| var_names start with `Zm` or `Zm00001` | Zea_mays |
| var_names start with `Os` or `LOC_Os` | Oryza_sativa |
| var_names start with `Brara\|Bra` | Brassica_rapa |
| var_names contain `Medtr\|MTR_` | Medicago_truncatula |
| var_names contain `Glyma\.` | Glycine_max |
| `adata.uns['organism']` set | Use that value directly |
| `adata.uns['species']` set | Use that value directly |

If detected, look up in `species-selection.md` and propose 1-3
model_species. **Always confirm with user** before running.

### Organ hint detection

| Signal | Detected organ |
|--------|----------------|
| h5ad path contains `leaf` | leaf |
| h5ad path contains `root` | root |
| `adata.uns['tissue']` or `adata.uns['organ']` set | Use that value |
| User prompt mentions tissue name | Match against `list_organs.py` output |
| Nothing detected | Default to the most common organ in the dataset, or ask |

### Cluster size for top_num estimation

Compute median cluster size `N_med`. Then:

```text
if N_med < 50:   top_num = 50
elif N_med < 200: top_num = 30   (default)
elif N_med < 1000: top_num = 20
else:            top_num = 15
```

Rationale: smaller clusters need more candidate genes to find
specific markers; larger clusters have abundant signal and can
afford smaller top_num.

## What to do with the inspection result

1. **Compose a config.yaml** with the inferred values (the script already
   recommends `cluster_key` and `top_num`; you decide `model_species` and
   `organ` from the species/organ hints).
2. **Present it to the user** with reasoning for each non-trivial choice
   (species, organ) — unless the user has asked you to run fully autonomously,
   in which case log the reasoning instead of blocking on confirmation.
3. **For technical parameters** (top_num, mode, decay_factor, etc.),
   just include them in the config with a brief inline comment.
4. **Do not re-read the h5ad.** Everything you need is in the
   `inspect_h5ad.py` output; a second `sc.read()` of a multi-GB file is the
   most common avoidable I/O cost (`mandatory-rules.md` §10).

Example presentation:

```
I inspected the h5ad file (one load, scripts/inspect_h5ad.py) and detected:
- cluster_key: "leiden" (15 clusters)
- species hint: var_names start with "Zm" -> Zea_mays (maize)
- organ: not detected in metadata
- recommended top_num: 30 (median cluster size = 180 cells)

Proposed config:
- model_species: [Zea_mays, Oryza_sativa, Sorghum_bicolor]
  (same family Poaceae, all 3 have gene models in xener)
- organ: leaf (default; please confirm)
- top_num: 30
```
