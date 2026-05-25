---
name: xener
description: |
  Xener is a cross-species single-cell cell type annotation tool using knowledge graphs.
  It maps marker genes from non-model species to model species via BLAST homology,
  then propagates weights through a knowledge graph to predict cell types per cluster.
  Supports sub-cluster refinement with Moran's I gene filtering.

  Use this skill when users want to annotate single-cell data, predict cell types,
  map genes across species, or run any xener-related tasks.

  Trigger whenever user mentions:
  - single-cell annotation, cell type prediction
  - xener, cross-species gene mapping
  - h5ad file cell type analysis
  - marker gene annotation for clusters
  - scRNA-seq cluster annotation

  CRITICAL: Do NOT write custom Python code. Use ONLY the provided scripts in `scripts/`.
  Do NOT import xener directly. Do NOT write your own pipeline code.
---

# Xener Skill

## MANDATORY RULES

**These rules are strict. Violations will cause incorrect behavior.**

1. **Use ONLY scripts** in `scripts/` directory to accomplish tasks
2. **Run scripts via bash/command line** — e.g., `python scripts/run_pipeline.py --config config.yaml`
3. **If a script fails, you MAY edit the script code** to fix the error
4. **Do NOT generate new code** outside of scripts. Only modify existing scripts when needed.
5. **NEVER import or use xener directly** — not in scripts, not in one-liners, not in any code

**This is absolutely forbidden**:
```bash
# NEVER do any of these:
python -c "from xener import Xener; ..."
python -c "import xener; ..."
from xener import Xener
import xener

# NEVER pass config via stdin:
python scripts/run_pipeline.py --config /dev/stdin
```

**Exception**: Only `python -c "import xener"` is allowed ONLY for checking if xener is installed. No other xener usage is permitted.

**Correct approach**:
```
python scripts/run_pipeline.py --config config.yaml
# If it fails, edit scripts/run_pipeline.py to fix it
python scripts/list_species.py
python scripts/list_organs.py
# If a script fails, edit the script to fix it
```

**Wrong approach**:
```bash
# Writing your own script to do the task
python my_custom_script.py
# NEVER write your own code
python -c "from xener import Xener; x = Xener(); x.list_species()"
```

## Installation

```bash
pip install xener
```

Or via the script:
```bash
python scripts/install.py
```

## Query Available Options

Reference species come from the BLAST database; organs come from the knowledge graph.

```bash
# List all available reference species (from BLAST database)
python scripts/list_species.py

# List all available organs (from knowledge graph)
python scripts/list_organs.py
```

**Always present results to the user and wait for confirmation** of the final `model_species` and `organ` before running the pipeline.

**Rule**: When showing species or organ lists to the user:
1. **Categorize** results (e.g., plants, animals, or by taxonomic group)
2. **Present by category** with a header for each group
3. **Format as 6-items-per-row aligned table** with centered content
4. **No truncation**, no "etc."

Example format:
```markdown
### Plants
| Species 1          | Species 2          | Species 3          | Species 4          | Species 5          | Species 6          |
|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|
| Arabidopsis_thaliana | Brassica_rapa      | Oryza_sativa       | Solanum_lycopersicum | Medicago_truncatula | Glycine_max        |

### Animals
| Species 1          | Species 2          | Species 3          | Species 4          | Species 5          | Species 6          |
|--------------------|--------------------|--------------------|--------------------|--------------------|--------------------|
| Homo_sapiens       | Mus_musculus       | Danio_rerio        |                    |                    |                    |
```

**Important**: You must categorize the results based on common biological knowledge. If you cannot determine the category, present them in a single group labeled "Other".

### Suggest Refinement Candidates

After step5, extract top cell types for agent to analyze and suggest refinement.

```bash
python scripts/suggest_refine.py \
    --input output/celltype_weight.zip \
    --topk 5 \
    --outdir output/
```

**Agent responsibilities**:

1. Run `suggest_refine.py` to get top-5 cell types per cluster
2. **Perform semantic deduplication**: identify cell types that refer to the same biological entity and group them (e.g., "phloem" + "vascular tissue" → same tissue type, keep only the higher-weight one)
3. **Check weight similarity**: for each cluster, compare the top 2 *distinct* cell types' `init_weight`. If `weight_2 / weight_1 > threshold` (default 0.5), the cluster is eligible for refinement
4. **Present results**: show eligible clusters with the deduplicated `candidate_celltypes` for user confirmation

**Example**:
```
Cluster 0:
  Top 5: phloem (0.8), vascular tissue (0.75), xylem (0.5), mesophyll (0.4), palisade (0.35)

  Agent deduplication:
  - phloem + vascular tissue = same vascular tissue → keep "phloem" (higher weight)
  - xylem = separate tissue
  - mesophyll + palisade = same tissue type → keep "mesophyll" (higher weight)

  Distinct types: phloem (0.8), xylem (0.5), mesophyll (0.4)
  Top 2: phloem (0.8) vs xylem (0.5), ratio=0.625 > 0.5 → ELIGIBLE
  Suggested: --cluster-id 0 --celltype phloem,xylem
```

## Available Scripts

All scripts are in the `scripts/` directory.

### Full Pipeline

```bash
python scripts/run_pipeline.py --config config.yaml
```

### Step-by-Step (for parameter tuning)

Each step saves intermediate results to `outdir/`. Rerun a step with different parameters to optimize.

```bash
# Step 1: Get marker genes
python scripts/step1_markers.py --input data.h5ad --cluster-key leiden --outdir output/

# Step 2: Calculate gene weights
python scripts/step2_weight.py --input output/marker_gene.zip --method prod --outdir output/

# Step 3: BLAST homology mapping
python scripts/step3_mapping.py --input output/marker_weight.zip \
    --fasta Arabidopsis_thaliana.fasta \
    --species Brassica_rapa \
    --weight-key pident \
    --pident 60 --evalue 0.05 --bitscore 200 \
    --mapping-strict 0 \
    --outdir output/

# Step 4: Get top-k genes
python scripts/step4_topk.py --input output/gene_homolo_weight.zip --k 30 \
    --multihomolo --outdir output/

# Step 5: Cell type annotation
python scripts/step5_annotate.py --input output/topk_markers.zip \
    --outdir output/annotation \
    --organ leaf \
    --mode path --decay-factor 0.7 \
    --ann-strict 0 \
    --candidate-annotation type1 type2
```

**Note**: `--candidate-annotation` restricts cell types. Run without it first to see available types.

### Step Parameter Reference

| Parameter | Step | Default | Description |
|-----------|------|---------|-------------|
| `--method` | step2 | `prod` | Weight calculation: `prod` or `sum` |
| `--weight-key` | step3 | `pident` | Homology weight column: `pident`, `evalue`, `bitscore` |
| `--pident` | step3 | `60` | Minimum percent identity filter |
| `--evalue` | step3 | `0.05` | Maximum e-value filter |
| `--bitscore` | step3 | `200` | Minimum bitscore filter |
| `--k` | step4 | `30` | Top N genes per cluster |
| `--multihomolo` | step4 | `true` | Keep multiple homologs per gene |
| `--mode` | step5 | `path` | Annotation mode: `node` (single type) or `path` (trajectory) |
| `--decay-factor` | step5 | `0.7` | Weight decay factor for graph propagation |
| `--organ` | step5 | `None` | Organ filter for knowledge graph |
| `--threshold` | step5 | `None` | Z-score threshold for cell type filtering |
| `--mapping-strict` | step3 | `0` | <0=ignore BLAST quality (all weights=1), 0=balanced, 1=suppress multi-copy gene families |
| `--ann-strict` | step5 | `0` | <0=exploratory (more cell types predicted), 0=balanced, 1=cleaner (1 type/marker), 2=strictest (1 type/cluster) |

### Sub-cluster Refinement

```bash
python scripts/refine_cluster.py \
    --input data.h5ad \
    --markers output/gene_homolo_weight.zip \
    --cluster-key leiden \
    --cluster-id 0 \
    --celltype type1,type2 \
    --organ leaf \
    --moran-i 0.5 \
    --split-method bindiv \
    --markergene-method diff \
    --strict 0 \
    --outdir output/
```

**Refinement parameters**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--moran-i` | `0.5` | Moran's I threshold for gene filtering. Range [-1, 1]. Closer to 1 = stricter. Invalid values skip filtering. |
| `--split-method` | `bindiv` | Cell assignment strategy: `bindiv` (binary division) or `argmax` |
| `--markergene-method` | `diff` | Marker gene set: `diff` (differential only) or `all` |
| `--key-added` | `xener_refine` | Column name in `adata.obs` for refined annotation |
| `--strict` | `0` | 0=default, >0=keep only max-confidence cell type per gene — cleaner sub-cluster annotation |

**Important**: `candidate_celltype` values must come from `celltype_weight.zip` for that specific cluster. Run step5 first to get valid cell type names.

## Config File Format

```yaml
# Required fields
cluster_key: leiden
model_species:
  - Brassica_rapa
non_model_fasta: Arabidopsis_thaliana.fasta
non_model_h5ad: data.h5ad
organ: leaf
outdir: output/

# Optional fields (defaults shown)
marker_weight_method: prod          # prod | sum
top_num: 30                         # Top N genes per cluster
homolo_weight_key: pident           # pident | evalue | bitscore
multihomolo: true                   # Keep multiple homologs per gene
decay_factor: 0.7                   # Weight decay for graph propagation
mode: path                          # node | path
threshold: null                     # Z-score threshold for filtering
mapping_strict: 0                   # <0=ignore BLAST quality, 0=balanced, 1=suppress multi-copy families
ann_strict: 0                       # <0=exploratory (more types), 0=balanced, 1=cleaner (1 type/marker), 2=strictest (1 type/cluster)
candidate_annotation:               # Restrict cell types
  - type1
  - type2
```

## Config File Validation Workflow

**Before running the pipeline, you MUST validate the config file:**

### Step 1: Find config file
First, check what config files exist in the current directory:
```bash
ls -la *.yaml *.yml 2>/dev/null | grep -v pip
```

**If multiple config files exist**: Ask user which one to use.
**If no config files exist**: Ask user whether they want to create a new one.
**If user provides a specific config path**: Use that path instead.

### Step 2: If config exists, validate it
1. **Parse the config file** and extract `model_species` and `organ` values
2. **Run validation scripts**:
   - `python scripts/list_species.py` — get list of valid species
   - `python scripts/list_organs.py` — get list of valid organs
3. **Check validity**:
   - **ALL** `model_species` values in config must appear in the species list (complete match, not partial)
   - The `organ` value must appear in the organ list (exact match)
   - All other required fields (`non_model_fasta`, `non_model_h5ad`, `outdir`) must be non-empty

### Step 3: If config is valid, proceed
If all validations pass, run:
```bash
python scripts/run_pipeline.py --config <config_file_path>
```

Replace `<config_file_path>` with the actual config file path (e.g., `config.yaml`, `my_config.yml`, etc.).

### Step 4: If config is missing or invalid
1. **If NOT_FOUND**: Ask user whether they want to create a new config file
2. **If INVALID** (bad species/organ/fields): Show the validation errors and ask user whether to fix the existing file
3. **If user confirms**: Create/fill config.yaml with placeholders:
   ```yaml
   cluster_key: leiden
   model_species: []
   non_model_fasta: ""
   non_model_h5ad: ""
   organ: ""
   outdir: output/
   candidate_annotation: []
   ```
4. Then query available options and guide user to fill in valid values
5. **Re-validate after changes** before running pipeline

**Never run `python scripts/run_pipeline.py --config /dev/stdin` or any stdin-based config passing.**

## Workflow

1. **Check if xener is installed** — run `pip show xener` or try `python -c "import xener"` (inline import only for checking)
2. **If NOT installed, install first** — run `pip install xener` or `python scripts/install.py`
3. **Validate config.yaml** using the steps above
4. Query available species (`list_species.py`) and organs (`list_organs.py`)
5. **Present options to user and wait for confirmation** of the final `model_species` and `organ` choices
6. Prepare config.yaml or run step-by-step
7. For parameter tuning: run individual steps with different parameters
8. After step5: run `suggest_refine.py` to get top-5 cell types per cluster
9. Agent performs semantic deduplication and checks weight similarity
10. Present refinement candidates to user for confirmation
11. Run `refine_cluster.py` for eligible clusters

**Important**: Always confirm species and organ with the user before proceeding. Wrong values will cause pipeline failure.

## Output Files

The pipeline generates the following outputs in `outdir/`:

| File | Description |
|------|-------------|
| `marker_gene.zip` | Marker genes per cluster |
| `marker_weight.zip` | Weighted marker genes |
| `blastp_{species}.zip` | BLASTP alignment results (one per model species) |
| `gene_homolo_weight.zip` | Homology-mapped genes with weights |
| `topk_markers.zip` | Top-k genes per cluster |
| `celltype_weight.zip` | All predicted cell types per cluster (use for refinement `candidate_celltype`) |
| `debug_params.yaml` | Actual parameters used in each step (for reproducibility) |
| `config.yaml` | Config copy (from `run_from_yaml` only) |
| `annotation/` | Per-cluster annotation XML files (`cluster_{id}_gene2celltype.xml`) |
| `refine_suggestions.json` | Suggested clusters for refinement (from `suggest_refine.py`) |

**Script intermediate files**: The step scripts save intermediate CSVs (e.g., `marker_gene.zip`, `marker_weight.zip`) for checkpointing. The full pipeline saves `.zip` files.

### debug_params.yaml

Records the actual parameter values used in each key step:

```yaml
cell_annotation:
  ann_strict: 0
  decay_factor: 0.7
  mode: path
  organ: leaf
  threshold: null
get_gene_weight:
  marker_weight_method: prod
get_topk_gene:
  multihomolo: true
  top_num: 30
mapping:
  bitscore: 200
  evalue: 0.05
  homolo_weight_key: pident
  mapping_strict: 0
  model_species:
  - Oryza_sativa
  pident: 60
```