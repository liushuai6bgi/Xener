# Refinement Workflow

After the main pipeline (Step 5), clusters may be "mixed" — the
top-1 predicted cell type has a weight close to the top-2.
Xener supports splitting such clusters into subtypes using
Moran's I gene filtering and binary / argmax cell assignment.

## When to suggest refinement

The agent (or user) should suggest refinement for a cluster if
**both** conditions hold after step 5:

1. **Semantic deduplication** identifies redundant top-5 cell
   types (e.g., "phloem" + "vascular tissue" → same tissue,
   keep the higher-weight one).
2. **Weight ratio** of top-2 distinct cell types:
   `weight_2 / weight_1 > 0.5` (configurable threshold).

## Step A: Get refinement candidates

```bash
python scripts/suggest_refine.py \
    --input output/celltype_weight.csv \
    --topk 5 \
    --outdir output/
```

→ writes `refine_suggestions.json`

The agent then performs:
- **Semantic dedup** — group synonyms (e.g., "phloem" + "vascular
  tissue" → keep "phloem")
- **Weight check** — for each cluster, compare top-2 distinct
  weights; if ratio > 0.5, mark as **eligible**

## Step B: Confirm with the user

Present eligible clusters with their deduplicated candidate cell
types and the suggested command. **Wait for confirmation** before
running refinement.

### Example agent reasoning

```
Cluster 0:
  Top 5: phloem (0.8), vascular tissue (0.75), xylem (0.5),
          mesophyll (0.4), palisade (0.35)

  Agent dedup:
  - phloem + vascular tissue → same → keep "phloem" (0.8)
  - mesophyll + palisade → same → keep "mesophyll" (0.4)

  Distinct top-2: phloem (0.8) vs xylem (0.5), ratio=0.625 > 0.5
  → ELIGIBLE

  Suggested: --cluster-id 0 --celltype phloem,xylem
```

## Step C: Run refinement

```bash
python scripts/refine_cluster.py \
    --input data.h5ad \
    --markers output/gene_homolo_weight.csv \
    --cluster-key leiden \
    --cluster-id 0 \
    --celltype phloem,xylem \
    --organ leaf \
    --moran-i 0.5 \
    --split-method bindiv \
    --markergene-method diff \
    --outdir output/
```

**Important**: `--celltype` values must come **exactly** from
`celltype_weight.csv` for the target cluster. Do not invent names.

## Step D: Inspect the refined annotation

Refinement writes a new column (default `xener_refine`) to
`adata.obs` containing the sub-cluster label. Visualize the
refined UMAP to confirm the split is biologically meaningful.

## Return values from `refine_single_cluster`

When you call the Python API (or inspect the script output),
`refine_single_cluster` returns **four** values:

| Return value | Type | Description |
|--------------|------|-------------|
| `geneCount` | list[(celltype, gene_count, gene_list)] | Markers per cell type |
| `diffgeneCount` | list[(celltype, diff_gene_count, diff_gene_list)] | Unique (differential) markers per cell type |
| `annotation` | pd.DataFrame | `adata_sub.obs[[key_added]]` — refined sub-cluster labels |
| `gene2celltype_g` | networkx.Graph | Gene → homolo → celltype graph for this refinement |

The CLI wrapper (`scripts/refine_cluster.py`) saves the first
three as a CSV and the fourth as a GEXF file:

```
output/refine_output/
├── refined_<cluster_id>.csv            # annotation DataFrame
└── refined_<cluster_id>_gene2celltype.gexf   # gene2celltype_g graph
```

The `gene2celltype_g` graph is the **provenance** of the
refinement: which marker genes drove the sub-cluster split,
through which homologs, into which cell types. Inspect it with
Gephi, Cytoscape, or `networkx`:

```python
import networkx as nx
g = nx.read_gexf('output/refine_output/refined_4_gene2celltype.gexf')
print(f"Nodes: {g.number_of_nodes()}, Edges: {g.number_of_edges()}")
# Inspect a specific gene's neighborhood:
neighbors = list(g.neighbors('AT1G52050'))
```

## Parameters

See `parameters.md` for the full refinement parameter table.
The three most important knobs:
- `--moran-i` (0.5) — gene spatial autocorrelation threshold;
  higher = stricter
- `--split-method` — `bindiv` (binary divide) vs `argmax`
  (assign by max score)
- `--markergene-method` — `diff` (differential only) vs
  `all` (all top-k genes)
