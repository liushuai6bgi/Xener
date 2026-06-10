# Mandatory Rules

These rules are strict. Violations will cause incorrect behavior or silent
errors. Read this file in full before any operation.

## 1. Use ONLY scripts in `scripts/`

**Why**: The Xener skill is a thin orchestration layer over a complex
bioinformatics toolchain. Bypassing scripts means losing the checkpointing,
parameter validation, and debug_params.yaml reproducibility record that the
official pipeline provides.

**Allowed**:
```bash
python scripts/run_pipeline.py --config config.yaml
python scripts/step1_markers.py --input data.h5ad --cluster-key leiden --outdir output/
python scripts/list_species.py
```

**Forbidden**:
```bash
python -c "from xener import Xener; ..."
python -c "import xener; ..."
python my_custom_script.py
```

**Exception**: Only `python -c "import xener"` is permitted to check
whether xener is installed. No other xener usage is allowed in one-liners.

## 2. Never pass config via stdin

```bash
# WRONG
python scripts/run_pipeline.py --config /dev/stdin < config.yaml

# CORRECT
python scripts/run_pipeline.py --config /path/to/config.yaml
```

**Why**: stdin-based config passing breaks reproducibility (no record of the
exact config used) and prevents `debug_params.yaml` from being written.

## 2b. Write `config.yaml` as UTF-8, and keep comments ASCII-safe

`run_pipeline.py` now reads the config as UTF-8 explicitly, so plain UTF-8
configs are safe. But be aware of the failure this guards against: on a
non-UTF-8 Windows console (e.g. zh-CN, default code page GBK), a tool that
opens the file with the *platform default* encoding will raise
`UnicodeDecodeError: 'gbk' codec can't decode byte ...` the instant a YAML
comment contains a non-ASCII character — the classic culprit is an **em-dash
(`—`)** or other "smart punctuation" pasted into a comment.

Practical guidance when composing a config:
- Prefer **ASCII** in comments (`-` instead of `—`, plain quotes). It is the
  one encoding that can never trip this.
- If you do write non-ASCII, ensure the file is saved as UTF-8.
- If you hit a `UnicodeDecodeError` mentioning `gbk`/`cp936`, the fix is the
  config file's bytes, not the pipeline — rewrite the offending comment as
  ASCII and re-run.

## 3. Always confirm `model_species` and `organ` with the user

Wrong values will cause pipeline failure. Before running, run:

```bash
python scripts/list_species.py
python scripts/list_organs.py
```

Present the available options to the user (categorized by Plants/Animals/etc.,
in 6-column tables — see main `SKILL.md` for the formatting rule), then
**wait for explicit confirmation** before running the pipeline.

## 4. Use script-based recovery, not custom code

If a script fails, you may edit the script to fix the bug, then re-run it.
**Do not write new scripts** to work around a failure.

## 5. Config field validation

Before running, the config must satisfy:
- `model_species`: list, each element must be in `scripts/list_species.py` output
- `organ`: string, must be in `scripts/list_organs.py` output
- `non_model_h5ad`, `non_model_fasta`, `outdir`: non-empty strings

See `workflows/config-validation.md` for the full validation procedure.

## 6. Refinement candidate celltypes must come from step 5 output

`--celltype` arguments to `scripts/refine_cluster.py` must be values that
appear in `celltype_weight.csv` for the target cluster. Do not invent cell
type names — they will fail at runtime.

## 7. Optional parameters inherit defaults

You do not need to specify every parameter. The pipeline reads defaults
from the package config. Override only when you have a reason.

## 8. Completeness over demonstration

When the user provides a real dataset and asks for cell-type
annotation, the task is **complete-annotation**, not a workflow
demo. The agent must refine **every** cluster where the top-2
distinct `init_weight` ratio is > 0.5 — not a representative
subset, not a "demonstrably mixed" subset, not a 1–3-cluster
cap. Refining N eligible clusters costs ~N×1 minute and produces
the per-cluster sub-population the user needs. The previous
"pick up to 3" heuristic was a demo-mode leftover and is
**withdrawn** — see `autonomous-decision-making.md` §"Completeness
vs. demonstration" for the full reasoning.

Specific anti-patterns that violate this rule:

- ❌ Selecting clusters "for lineage coverage" (one vasculature,
  one ground tissue, one root cap) — that is a presentation
  choice, not an annotation choice.
- ❌ Skipping a high-ratio cluster because the top-2 cell types
  "look biologically unrelated" (e.g. quiescent center + root
  hair cell) — that unrelatedness is the *reason* to refine; it
  is a signal that the clustering may be bad or that KG
  propagation biased step 5. The refinement is the test.
- ❌ Treating 3/3 uniform refinements as "sufficient evidence" the
  clustering is good — that is a sample of 3, not the dataset.
- ❌ Capping at any small number to "leave time for visualization"
  or "keep the narrative clean". The user does not pay for agent
  time; do the work.
