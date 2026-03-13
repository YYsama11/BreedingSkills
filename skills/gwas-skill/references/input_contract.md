# Input Contract

## Required files

Default layout assumes a working directory with a `data/` folder containing:

- `sample_ids.txt` — target sample IDs, one per line
- `phenotype_matrix.tsv` — phenotype table with metadata columns followed by replicate columns such as `SampleA_1`
- `genotype_panel.tped` and `genotype_panel.tfam` — genotype prefix
- `kinship_matrix.tsv` — kinship matrix matching the full `.tfam` order
- `covariates.tsv` — covariates matching the full `.tfam` order
- `annotation.gff3` — annotation GFF
- `reference.fa.fai` — reference index for chromosome lengths

## Important assumptions

- Replicates are identified by suffixes `_1`, `_2`, `_3`.
- Replicates do not need to be adjacent; grouping is done by sample ID prefix.
- `EMMAX` uses `tped/tfam`, kinship, and phenotype files aligned to the subset `.tfam` order.
- The candidate-gene workflow expects annotation coordinates in the same reference system as the genotype SNP positions.

## Environment variables

Set these when filenames or tool paths differ from defaults:

- `GWAS_WORKSPACE_DIR` — root analysis directory
- `GWAS_DATA_DIR` — input data directory, default `${GWAS_WORKSPACE_DIR}/data`
- `GWAS_ID_FILE`
- `GWAS_LIPID_FILE`
- `GWAS_GENO_PREFIX`
- `GWAS_KINSHIP_FILE`
- `GWAS_COVARIATE_FILE`
- `GWAS_GFF_FILE`
- `GWAS_FASTA_FAI_FILE`
- `GWAS_THREADS`
- `GWAS_EMMAX_PARALLEL`
- `PLINK_BIN`
- `EMMAX_BIN`

## Workspace model

- The skill repository contains the scripts.
- `GWAS_WORKSPACE_DIR` points to the project workspace that contains `data/` and receives `analysis/`.
- You should run `bash scripts/run_full_pipeline.sh` from inside the `skills/gwas-skill/` directory, while `GWAS_WORKSPACE_DIR` points to the external workspace.

## Minimal launch example

```bash
export GWAS_WORKSPACE_DIR=/path/to/project
export GWAS_DATA_DIR=/path/to/project/data
export GWAS_THREADS=6
export GWAS_EMMAX_PARALLEL=10
bash scripts/run_full_pipeline.sh
```
