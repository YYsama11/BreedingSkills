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
- `analysis/gwas/emmax/results/*.log`

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

## What this workflow deliberately excludes

The core GWAS skill does not include:

- raw-read alignment
- variant calling
- SNP functional annotation
- LD/QTL analysis
- candidate-gene extraction
- phenotype-domain-specific subclassification

Those belong in separate skills.
