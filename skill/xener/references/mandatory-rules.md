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

**Note on `scripts/_xener_init.py`**: this is a *helper module*, not a CLI —
the other scripts import `build_xener()` from it so they all construct the
`Xener` object the same way (and can honor an optional `--init-config`). Do not
run it directly, and do not hand-write your own `Xener(...)` construction to
point at a custom KG / BLAST database. The supported way to use custom
infrastructure is the `--init-config` flag plus `scripts/init_xener.py` to
validate it — see `workflows/initialization.md`.

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

## 9. Long-running steps: wait by notification, NEVER by `sleep`

The two long steps are the **pipeline** (`run_pipeline.py`) and **batch
refinement** (`refine_cluster.py --plan`). With a warm BLAST cache they finish
in a few minutes; plotting is faster still. How you wait for them matters —
this rule exists because an earlier run wasted ~11 minutes polling with
`sleep N && grep`.

**Run them one of two ways:**

1. **Foreground with a generous Bash `timeout`** (e.g. `timeout: 600000` = the
   10-minute max). Simplest; the call returns when the step is done.
2. **Background** (`run_in_background: true`). The harness **re-invokes you
   automatically when the task completes** — so after launching, *end your
   turn* and do nothing until that completion notification arrives.

**NEVER do this:**

```bash
# FORBIDDEN: blind sleep-polling
sleep 90 && grep ... outdir/xener.log
```

Why it is wrong:
- `sleep` **blocks your turn**. The completion notification can only be
  processed when you are idle, so sleeping *delays your own reaction* to the
  task finishing.
- It is a **blind timer**: you either sleep too short (read a half-written
  log, then sleep again — wasted round-trips) or too long (the task finished
  long ago and you are still asleep). A timer cannot align with an event.
- `sleep ≥ 120` additionally **collides with the Bash command timeout**
  (default 120 000 ms) — the call is killed with exit 143 and you have waited
  two minutes for nothing. This actually happened.

**Do this instead:**

| Goal | Correct mechanism |
|------|-------------------|
| Just wait for completion | Background the task, then **end your turn**. The completion notification re-invokes you. |
| Wait synchronously, in-turn | `TaskOutput` with `block: true` — returns the instant the task ends, not on a guessed timer. |
| Peek at progress mid-run | A **bare** `grep outdir/xener.log` with **no** `sleep` prefix. It is instant; an empty result just means "nothing logged yet", which is itself information. Do not wrap it in a poll loop. |
| Come back later without blocking | `ScheduleWakeup` (non-blocking timer), then **end your turn**. Unlike `sleep`, it does not occupy you and does not mask the completion notification. |

The unifying principle: **status checks are instantaneous and need no timer;
waiting is event-driven (notification / `TaskOutput` / `ScheduleWakeup`), never
`sleep`-driven.**

## 10. Inspect the h5ad ONCE; never re-read it when a CSV will do

The target `.h5ad` is frequently multi-GB and a single read (~30-45 s for ~1 GB)
is the **dominant I/O cost of the whole workflow**. Treat each read as expensive.

- **Inspection is a single command:** `python scripts/inspect_h5ad.py <h5ad>`.
  It loads the file once and prints everything needed to compose a config
  (cluster_key, species hint, organ hint, cluster sizes, recommended `top_num`,
  `.X`/raw sanity). **Do not** follow it with ad-hoc `python -c "import scanpy;
  sc.read(...)"` reads to fetch "one more field" — every field is already in
  its output. (This is the *only* sanctioned inline-scanpy path besides rule §1.)
- **The quality gate does not re-read the h5ad.** `run_pipeline.py` calls the
  gate in-process and the gate derives cluster sizes from the lightweight
  `{dataset}_annotation.csv`, not the h5ad.
- **Plotting does not re-read the h5ad.** `plot_umap.py` takes the same
  annotation CSV (per-cell UMAP coords + labels).

If you find yourself about to call `sc.read()` on the big file a second time,
stop: the answer is almost certainly already in `inspect_h5ad.py`'s output or
in a `*_annotation.csv` / `*.csv` already on disk.

## 11. Hard constraints outrank the quality gate

**A value fixed by the biology of the dataset, by the user, or by the
available infrastructure is a HARD CONSTRAINT. The post-run quality gate
(`check_output.py`, `self-tuning-protocol.md`) is a SOFT SIGNAL. Never change
a hard constraint to make a soft signal pass. If the only way to pass the gate
is to alter a hard constraint, the correct outcome is to ACCEPT the gate
failure and document why — not to alter the constraint.**

This rule is data-independent; it applies to every run.

### What is hard vs soft

| HARD constraint — do NOT change to pass the gate | Why it is hard |
|---|---|
| `organ`, once it reflects a known/confirmed tissue of the sample | Biological fact about the data |
| `cluster_key` | What the h5ad actually contains |
| input `non_model_h5ad` / `non_model_fasta` | The data under analysis |
| `model_species` **phylogenetic boundary** (must be genuine relatives of the target) | Validity of homology transfer |
| any value the **user explicitly set or confirmed** | User intent outranks the agent |
| `model_species` / `organ` membership in the KG + BLAST DB | Infrastructure fact; an absent value cannot run |

| SOFT parameter — tune these freely to improve the gate |
|---|
| `model_species` membership *inside* the valid boundary (how many relatives; whether to add a KG-rich close relative) |
| BLAST thresholds: `pident`, `evalue`, `bitscore` |
| `top_num`, `multihomolo`, `homolo_weight_key`, `marker_weight_method` |
| `mode`, `decay_factor`, `threshold`, `ann_strict`, `mapping_strict` |
| `candidate_annotation` (restrict the output cell-type set) |

The split on `model_species` is the subtle part: *which clade* is hard (biology);
*how many members within that clade* is soft (tuning). Adding a closer relative
that already exists in the KG is a soft fix. Switching to a phylogenetically
wrong species, or turning the organ filter off, is changing a hard constraint.

### Protocol when the gate fails

1. Identify which check failed (`self-tuning-protocol.md`) and the exact cause
   (`log-interpretation.md`).
2. Ask: **can a SOFT parameter fix it?** (add a valid closer relative, relax
   BLAST thresholds, raise `top_num`, drop `threshold`, restrict
   `candidate_annotation`, ...). If yes, do that and re-run.
3. If the only fix would change a HARD constraint — e.g. "switch `organ` away
   from the sample's real tissue", "set `organ=None`/`Unknown` purely to drop
   the filter", "add a phylogenetically unrelated species just for coverage" —
   **STOP. Do not do it.** Keep the hard constraint, accept the gate failure,
   and write the justification to `outdir/autonomous_log.md`: which check
   failed, which soft fixes you tried and their effect, and why no soft fix
   exists.
4. A gate failure that survives step 3 is an **accepted failure, not a failed
   run.** The deliverable is the biologically-correct annotation plus the
   documented residual gate signal — NOT a gate-green run built on a falsified
   constraint.

### The organ trap (most common instance)

When check 1 (KG miss) fails, the gate message and some docs say
"model_species too narrow for the chosen organ." Only these responses are legal:

- ✅ Add a closer relative that exists in the KG (soft), OR relax BLAST, OR —
  **only if `organ` was an unconfirmed guess** — reconsider the organ.
- ❌ Change `organ` to a tissue the sample is NOT, or to `None`/`Unknown`
  purely to lift the filter, **when the organ is a confirmed property of the
  data.** A high KG miss under the *correct* organ usually means the KG is
  simply shallow for that (clade x organ) pair — a property of the KG, not a
  config error. The honest result is an accepted gate failure under the correct
  organ.

**Diagnostic tell:** if you applied the gate-suggested soft lever (added a
KG-rich relative) and the metric barely moved, the bottleneck is NOT the soft
lever — it is a hard constraint's coverage in the KG. That is your cue to
ACCEPT, not to start changing hard constraints.

### Red flags — STOP, you are about to violate this rule

- "organ=None / organ=Unknown is THE FIX — it directly attacks the KG-miss check."
- "The docs say 'try a different organ', so changing the tissue is sanctioned."
- "A failing gate means the run isn't done; I must make it pass."
- "I correctly diagnosed that organ is the binding lever, therefore I should change organ."
- "Root passes the gate, so I'll annotate this stem sample against Root."

All of these mean: you are trying to pass a SOFT signal by breaking a HARD
constraint. Keep the constraint. Accept and document the gate result.
