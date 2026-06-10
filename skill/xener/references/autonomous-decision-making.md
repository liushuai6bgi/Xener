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
| Refinement targets | ✅ Yes, run on **all** eligible | After step 5, run `scripts/suggest_refine.py --topk 5` to get top-5 per cluster. Apply semantic dedup (e.g. "phloem" + "vascular tissue" → "phloem"). Mark cluster as eligible if top-2 distinct init_weight ratio > 0.5. **Refine every eligible cluster, not a subset.** The previous "pick up to 3" guidance was a demo-time heuristic and produced incomplete annotations — see `mandatory-rules.md` §8. Prioritization is *only* for ordering (tightest ratio first, then largest cluster first), never for *excluding* clusters. Do not pre-filter by "biological lineages with known subtypes" — refine cross-lineage high-ratio clusters (e.g. quiescent center + root hair cell at ratio 1.0) too; that pairing is a signal the clustering may be bad, not a reason to skip. Run `scripts/refine_cluster.py` on every eligible cluster; do not ask the user. |
| Refinement method (`--markergene-method`, `--split-method`) | ✅ Yes | **Use `--markergene-method all --split-method argmax`** — empirically the highest refinement success rate. The script defaults (`diff` / `bindiv`) frequently fail to split: `diff` returns 0 markers when the two candidates share KG homologs, and `bindiv` then assigns all cells to one candidate (`[N, 0]`). Only fall back to defaults if you specifically need stricter differential/binary behavior. |
| Refinement `--markers` input | ✅ Yes | **Always pass `gene_homolo_weight.csv` (Step 3), never `topk_markers.csv` (Step 4).** Refinement collects the cluster's homolog set from this file to query the KG; `topk_markers.csv` is truncated to `top_num` genes/cluster, which starves the split and makes nearly every cluster refine to one uniform label (a silent, exit-0 failure that *looks like* clean clustering). Only `group`/`gene`/`homolo` are read — coverage matters, weights do not. |

## Decision protocol for the LLM

1. **Read user prompt carefully.** Extract any explicit values for
   `organ`, `model_species`, `cluster_key`, etc. If the user provides
   a real dataset and asks for "annotation" (without further
   qualification), treat the task as **complete-annotation** mode —
   see "Completeness vs. demonstration" below — and plan to refine
   every eligible cluster, not a representative subset.
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

## Completeness vs. demonstration

The skill serves two distinct user intents. The agent must
distinguish them at the start of the task and act accordingly:

- **Complete-annotation mode** (default for any real dataset the
  user hands over): the user wants the dataset fully annotated.
  Refine every eligible cluster, integrate results back into the
  h5ad, and present a single, merged `xener_refine` column. The
  artifact of "I refined clusters 6, 15, 27 and they were uniform"
  is not complete — the user needs to know what happened to all
  33 eligible clusters. Refinement of N clusters costs ~N×1 min,
  which is acceptable work, not a budget to manage.
- **Demonstration mode** (only when the user explicitly asks "show
  me how this works", "demo", or "walk me through one example"):
  the user wants to see the workflow. Then it's fine to refine
  1–3 representative clusters to keep the run short.

**If unclear, default to complete-annotation mode.** "Annotation
of this dataset" is not a demo request.

## Anti-patterns to avoid

- ❌ **Asking the user about every parameter.** This defeats the
  purpose of an autonomous agent. Use the heuristics above.
- ❌ **Hard-coding default values for biological parameters.** Species
  and organ are dataset-specific. Always reason.
- ❌ **Running without inspecting the h5ad first.** You'll waste time
  on the wrong cluster_key or species.
- ❌ **Re-running the entire pipeline when only step 5 needs tweaking.**
  Use the checkpoint system — re-run from step 3 / 4 / 5 only.
- ❌ **Refining only a "representative" subset of eligible clusters.**
  This is a demo-mode leftover and produces incomplete annotations.
  Refine every eligible cluster. See "Completeness vs. demonstration".
- ❌ **Skipping cross-lineage high-ratio clusters on the grounds that
  the top-2 cell types "look biologically unrelated".** That
  unrelatedness is the *reason* to refine — it is a signal that
  the clustering may be wrong, the cluster may be a doublet, or
  KG propagation may have biased step 5. Refinement is the test,
  not a thing you do only after the test passes.
- ❌ **Picking "representative" lineage coverage (e.g. "one vasculature
  + one root cap + one ground tissue") to reduce the number of
  refinements.** Coverage is a presentation choice, not an
  annotation choice. Every eligible cluster gets refined.

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
