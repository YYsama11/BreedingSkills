# Workflow

## Phase 1 — Prepare trait-wise EMMAX phenotype files

Script:

- `scripts/prepare_emmax_inputs.py`

Inputs:

- trait matrix
- optional trait metadata
- covariate table

Outputs:

- `analysis/emmax_inputs/trait_manifest.tsv`
- `analysis/emmax_inputs/phenotypes/*.tsv`
- `analysis/emmax_inputs/covariates_emmax.tsv`

## Phase 2 — Run EMMAX across traits

Scripts:

- `scripts/run_emmax_batch.py`
- `scripts/run_emmax_trait.sh`

Inputs:

- trait manifest
- `tped/tfam`
- kinship matrix
- covariates

Outputs:

- `analysis/gwas/emmax/results/*.ps`
- `analysis/gwas/emmax/results/*.reml`

## Phase 3 — Build SNP reference arrays

Script:

- `scripts/prepare_snp_reference.py`

Inputs:

- `TPED` or `BIM`
- optional `FAI`

Outputs:

- `analysis/gwas/reference/chrom.npy`
- `analysis/gwas/reference/pos.npy`
- `analysis/gwas/reference/cum_pos.npy`
- `analysis/gwas/reference/chrom_ticks.tsv`

## Phase 4 — Summarize GWAS results

Script:

- `scripts/summarize_emmax_results.py`

Outputs:

- Manhattan plots
- QQ plots
- per-trait GWAS summary table
- significant SNP table
- shared SNP summary
- pleiotropic SNP table

## Optional standalone plotting

- `scripts/plot_manhattan_qq.R`
- `scripts/run_gwas_plot.sh`

This optional layer renders Manhattan + QQ directly from a generic GWAS result table.

## Logging behavior

- `analysis/logs/runtime/command_manifest.tsv` is always kept
- per-command runtime `*.cmd.log` and `*.top.log` are only kept if `KEEP_RUNTIME_LOGS=true`
- tool-generated `*.log` files are removed after successful execution unless `KEEP_TOOL_LOGS=true`

## What this workflow deliberately excludes

The core GWAS skill does not include:

- raw-read alignment
- variant calling
- SNP functional annotation
- LD/QTL analysis
- candidate-gene extraction
- phenotype-domain-specific subclassification

Those belong in separate skills.
