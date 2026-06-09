# Config Validation Workflow

Before running the pipeline, you MUST validate the YAML config.
This prevents runtime failures from bad species/organ names or
missing fields.

## Step 1: Find the config file

```bash
ls -la *.yaml *.yml 2>/dev/null | grep -v pip
```

- **Multiple configs found**: ask the user which one to use.
- **No config found**: ask whether to create one. The template is
  in `config-schema.md` (Minimal example).
- **User specifies a path**: use that path directly.

## Step 2: Parse the config

Read the file and extract `model_species` (list) and `organ`
(string). The skill should also check that the other required
fields are non-empty:
- `non_model_fasta`
- `non_model_h5ad`
- `outdir`
- `cluster_key`

## Step 3: Validate against available options

```bash
python scripts/list_species.py
python scripts/list_organs.py
```

**Check `model_species`**:
- Every element of `model_species` in the config must appear
  **exactly** (case-sensitive) in the species list.

**Check `organ`**:
- The `organ` value must appear **exactly** in the organ list.

**Present the categorized lists** to the user (see the formatting
rule in `SKILL.md` — 6-column tables grouped by Plants / Animals /
etc.), then **wait for explicit confirmation** before running.

## Step 4: Proceed or fix

- **All checks pass** → run `scripts/run_pipeline.py --config <path>`
- **`model_species` invalid** → show the user the valid options,
  ask which to use
- **`organ` invalid** → show the user the valid options, ask
  which to use
- **Required field empty** → ask the user to fill it in

## Step 5: After fixing, re-validate

If the user edits the config, re-run the full validation
(`list_species.py` + `list_organs.py` + field checks) before
running the pipeline. Do not skip re-validation.

## Forbidden

```bash
# NEVER pass config via stdin
python scripts/run_pipeline.py --config /dev/stdin
```
