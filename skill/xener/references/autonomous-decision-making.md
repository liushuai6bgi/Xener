# Autonomous Decision-Making Guide

This file defines **when and how the LLM agent should make decisions on its
own** versus when it should ask the user. The goal is to maximize the
agent's autonomy on technical parameters while keeping the user in control
of biological / experimental choices.

## Decision matrix

| Parameter | Agent decides alone? | Agent's heuristic |
|-----------|----------------------|-------------------|
| `cluster_key` | ✅ Yes | Scan `adata.obs.columns` for clustering results. Prefer `leiden` → `louvain` → any column with low cardinality. |
| `organ` | ✅ Yes, with reasoning | Infer from h5ad metadata (`adata.uns`, var names) or user prompt keywords. Always confirm with user before running. |
| `model_species` | ✅ Yes, with reasoning | If h5ad species is identifiable: pick 1-3 phylogenetically close species from `list_species.py` (see species-selection.md). **If the target species is itself a well-annotated model organism in `list_species.py` (e.g., Arabidopsis, human, mouse, zebrafish), include it in `model_species` even though the config field is `non_model_*` — the field name is a UX convention, not a biological constraint.** Self-mapping maximizes KG coverage; the cross-species independence concern is secondary and can be addressed by adding 1-2 close relatives on top. Else: ask the user. |
| `marker_weight_method` | ✅ Yes | Default `prod`. Switch to `sum` if log fold-change distribution is heavily right-skewed. |
| `top_num` | ✅ Yes, dynamic | `max(20, min(50, n_marker_genes_per_cluster // 2))`. Re-run with higher value if Step 5 returns empty. |
| `homolo_weight_key` | ✅ Yes | Default `pident`. Switch to `bitscore` if pident distribution has a long tail (no clear cutoff). |
| `multihomolo` | ✅ Yes | `true` for cross-species (typical), `false` for very close species with 1:1 orthology. |
| `pident` | ✅ Yes, iterative | Start at 60, lower to 40 if too few homologs. Raise to 75 for closely related species. |
| `evalue` | ✅ Yes | Default 0.05. Raise to 1e-3 only if homolog count is very low. |
| `bitscore` | ✅ Yes | Default 200. Lower to 100 if homolog count is very low. |
| `threshold` | ✅ Yes, data-driven | Compute z-score of `init_weight` across cell types in `celltype_weight.csv`. Set threshold = 1.5 if distribution is clean, else null. |
| `mode` | ✅ Yes | Default `path`. Switch to `node` if top-1 weight is >5x top-2 (clean annotation). |
| `decay_factor` | ✅ Yes | Default 0.7. Lower (0.5) for sparser graphs, raise (0.9) for dense ones. |
| `mapping_strict` | ✅ Yes | Default 0. Set to 1 if multiple model_species produce conflicting top homologs. Set to -1 (loose) for novel species. |
| `ann_strict` | ✅ Yes | Default 0. Set to 1 if multiple clusters share identical top cell types. |
| `num_threads` | ✅ Yes | Read `os.cpu_count()`. Use max(1, n_cpu - 2). |
| `outdir` | ⚠️ Ask user if not obvious | Default: `./xener_output/<dataset_basename>`. |
| Refinement targets | ✅ Yes, with reasoning | After step 5, run `scripts/suggest_refine.py --topk 5` to get top-5 per cluster. Apply semantic dedup (e.g. "phloem" + "vascular tissue" → "phloem"). Mark cluster as eligible if top-2 distinct init_weight ratio > 0.5. Pick up to 3 eligible clusters, prioritizing (a) tight ratios and (b) biological lineages with known subtypes (epidermis, ground tissue, vasculature, immune). Run `scripts/refine_cluster.py` directly — do not ask the user. |

## Decision protocol for the LLM

1. **Read user prompt carefully.** Extract any explicit values for
   `organ`, `model_species`, `cluster_key`, etc.
2. **Inspect the h5ad file** with `python -c "import scanpy as sc; adata =
   sc.read('<path>'); print(adata.obs.columns); print(adata.uns.keys())"`.
3. **Decide parameters** using the heuristics above. Write the
   decision + reasoning in your scratchpad.
4. **For biological parameters** (organ, model_species): still ask
   the user to confirm, but **propose 1-3 candidates** with reasoning
   instead of asking open-ended. The user should say "yes" or
   "use X instead" — not have to research from scratch.
5. **For technical parameters** (top_num, pident, mode, etc.): decide
   alone, and **log the decision** in the agent's scratchpad so the
   user can review.
6. **After running**, if any cluster's annotation looks suspicious
   (e.g., top-1 weight < 0.3, or many "unknown" annotations), the
   agent should **automatically re-run with adjusted parameters** —
   see `self-tuning-protocol.md` (workflows/).

## Anti-patterns to avoid

- ❌ **Asking the user about every parameter.** This defeats the
  purpose of an autonomous agent. Use the heuristics above.
- ❌ **Hard-coding default values for biological parameters.** Species
  and organ are dataset-specific. Always reason.
- ❌ **Running without inspecting the h5ad first.** You'll waste time
  on the wrong cluster_key or species.
- ❌ **Re-running the entire pipeline when only step 5 needs tweaking.**
  Use the checkpoint system — re-run from step 3 / 4 / 5 only.

## Self-tuning loop (optional but recommended)

After the initial run, the agent can:
1. Read `celltype_weight.csv` and check distribution of `init_weight`.
2. If median < 0.3 → annotations are weak. Re-run step 5 with
   `--threshold null --decay-factor 0.5` to widen the search.
3. If multiple clusters share the same top cell type → re-run step 5
   with `--ann-strict 1` to disambiguate.
4. If `gene_homolo_weight.csv` has very few rows → relax BLAST
   thresholds and re-run from step 3.

Document each re-run decision in `outdir/autonomous_log.md` for
reproducibility.
