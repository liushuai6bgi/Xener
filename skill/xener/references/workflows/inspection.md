# h5ad Inspection Workflow

Before running the pipeline, the agent **must inspect the h5ad file**
to extract metadata and make informed parameter decisions. This file
documents the inspection procedure.

## Why inspect first?

- The agent needs `cluster_key` — usually stored in `adata.obs`
- The agent may need species hints (var name prefix, uns metadata)
- The agent may need organ hints (uns metadata, file path)
- The agent needs to estimate `top_num` based on cluster sizes

## Inspection script (the ONLY inline import allowed)

```python
import scanpy as sc
adata = sc.read('<path_to_h5ad>')

# 1. Cluster columns
print("=== obs columns ===")
print(adata.obs.columns.tolist())
print("=== cardinality ===")
for col in adata.obs.columns:
    n_unique = adata.obs[col].nunique()
    if 2 <= n_unique <= 100:  # likely a clustering result
        print(f"  {col}: {n_unique} unique values")

# 2. Metadata
print("=== uns keys ===")
print(list(adata.uns.keys()))

# 3. Species hints from var_names prefix
if len(adata.var_names) > 0:
    sample = adata.var_names[:20].tolist()
    print("=== var_names sample ===")
    print(sample)

# 4. Basic shape
print(f"Shape: {adata.shape}")
print(f"obs: {adata.n_obs} cells, var: {adata.n_vars} genes")
```

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

1. **Compose a config.yaml** with the inferred values.
2. **Present it to the user** with reasoning for each non-trivial
   choice (species, organ).
3. **For technical parameters** (top_num, mode, decay_factor, etc.),
   just include them in the config with a brief inline comment.
4. **Wait for user confirmation** before running.

Example presentation:

```
I inspected the h5ad file and detected:
- cluster_key: "leiden" (15 clusters)
- species hint: var_names start with "Zm" → Zea_mays (maize)
- organ: not detected in metadata

Proposed config:
- model_species: [Zea_mays, Oryza_sativa, Sorghum_bicolor]
  (same family Poaceae, all 3 have gene models in xener)
- organ: leaf (default; please confirm)
- top_num: 30 (median cluster size = 180 cells)

Proceed? (y/n)
```
