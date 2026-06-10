# Refinement Workflow

After the main pipeline (Step 5), clusters may be "mixed" — the
top-1 predicted cell type has a weight close to the top-2.
Xener supports splitting such clusters into subtypes using
Moran's I gene filtering and binary / argmax cell assignment.

## When to suggest refinement

The agent (or user) should refine a cluster if the following
condition holds after step 5:

1. **Weight ratio** of top-2 distinct cell types (after semantic
   dedup — substring/containment, e.g. "phloem" vs "vascular
   tissue", or biological-synonym, e.g. "sieve element" vs
   "phloem") is `weight_2 / weight_1 > 0.5` (configurable
   threshold).

A high ratio (close to 1.0) means the top-2 cell types are
nearly tied in KG support — either a real sub-population is
hidden in the cluster, or step 5's KG propagation created
spurious secondary weight. **Either way, refine the cluster**
and let the result decide.

**Do not pre-filter eligible clusters by "biological plausibility
of the top-2 pair".** A cross-lineage pair (e.g. quiescent center
+ root hair cell at ratio 1.0) is *more* reason to refine, not
less — it is a signal that the clustering may be wrong, that
the cluster is a doublet, or that KG propagation biased step 5.
The agent's job here is to test, not to assume.

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

## Step B: Run refinement on EVERY eligible cluster

**In manual mode** — present the eligible list with their
deduplicated candidate cell types and the suggested commands.
**Wait for confirmation** before running refinement.

**In autonomous mode (complete-annotation)** — run refinement
on every eligible cluster. Order the runs by descending ratio
(tightest first), then descending cluster size. Do not skip
clusters, do not pick representatives, do not cap at 3. The
cap was a demo-time heuristic and is withdrawn — see
`mandatory-rules.md` §8.

**In autonomous mode (demonstration only — when the user
explicitly asks "show me how this works")** — it is acceptable
to run on 1–3 representative clusters. This is the *only*
context where "up to 3" applies.

Log the decision (which mode you are in, which clusters were
refined, in what order) in `outdir/autonomous_log.md` alongside
any other autonomous choices.

### Example agent reasoning (autonomous, complete-annotation)

```
33 / 36 clusters have top-2 / top-1 init_weight ratio > 0.5:
  - c6, c15, c27 (ratio ≥ 0.95) — known biological subtype pairs
  - c3, c20 (root cap subtypes) — same lineage as c15
  - c14, c34 (ratio = 1.00) — cross-lineage pairs, refine to
    test if clustering is bad
  - ... 25 more

Plan: refine all 33. Order: tightest ratio first.
Cluster c3 (lateral root cap + columella root cap cell, ratio 1.00)
  --celltype lateral root cap,columella root cap cell
Cluster c14 (quiescent center + root hair cell, ratio 1.00)
  --celltype "quiescent center","root hair cell"
...
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
    --split-method argmax \
    --markergene-method all \
    --outdir output/
```

**Important**: `--celltype` values must come **exactly** from
`celltype_weight.csv` for the target cluster. Do not invent names.

> **`--markers` MUST be `gene_homolo_weight.csv` (Step 3 output), NOT
> `topk_markers.csv` (Step 4 output).** This is the single most common
> refinement mistake and it fails *silently*: the run completes exit-0 but
> the cluster does not split (every cell gets one label).
>
> Why: `refine_single_cluster` uses the `--markers` table only to collect
> the cluster's **homolog set**, which it feeds to the KG to find markers for
> each candidate cell type. `topk_markers.csv` is truncated to `top_num`
> genes per cluster (e.g. 20), so the homolog set is thin, the KG returns too
> few candidate-discriminating markers, and `argmax` dumps all cells into one
> candidate → a `[N, 0]` result and a single-label cluster. `gene_homolo_weight.csv`
> carries the full BLAST homolog set per cluster, which is what the split needs.
> (The weight columns are ignored either way — only `group`/`gene`/`homolo`
> are read — so this is purely about homolog *coverage*.)
>
> Symptom to watch for: nearly every "eligible" cluster refines to a single
> uniform subtype. That is usually *this bug*, not evidence that the
> clustering is clean. `refine_cluster.py` prints a `[WARN] ... only N unique
> homologs` line when the `--markers` file looks truncated.

> **Use `--markergene-method all --split-method argmax`.** This combination
> has the **highest refinement success rate** in testing — it actually splits
> mixed clusters into sub-populations. The script defaults
> (`--markergene-method diff --split-method bindiv`) frequently produce *no*
> split: `diff` returns 0 differential markers whenever the two candidate
> cell types map to overlapping homologs in the KG (very common for close
> sub-types such as `lateral root cap` vs `columella root cap cell`), and
> `bindiv` then dumps every cell into one candidate, yielding a `[N, 0]`
> result and a single-label cluster. If you see `Diff gene counts: []` or
> `refine_single_cluster X total: ['<one type>']` for clusters you expected
> to split, switch to `all` + `argmax` and re-run.

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
- `--markergene-method` — `all` (all top-k genes) vs `diff` (differential
  only). **Use `all`**: `diff` returns 0 markers when candidates share KG
  homologs, which blocks the split.
- `--split-method` — `argmax` (assign by max score) vs `bindiv` (binary
  divide). **Use `argmax`**: it splits reliably; `bindiv` tends to collapse
  to one candidate.
- `--moran-i` (0.5) — gene spatial autocorrelation threshold;
  higher = stricter.

The empirically best combination for splitting mixed clusters is
**`--markergene-method all --split-method argmax`** (highest success rate).
