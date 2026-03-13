# Input Contract

## Required inputs

The core GWAS skill assumes the following required files:

- `trait_matrix.tsv`
- `genotype_panel.tped`
- `genotype_panel.tfam`
- `kinship_matrix.tsv`
- `covariates.tsv`

## Optional inputs

- `trait_metadata.tsv`
- `reference.fa.fai`

## Trait matrix format

The phenotype matrix must contain:

- column 1: `FID`
- column 2: `IID`
- column 3 onward: one or more quantitative traits

Example:

```text
FID   IID   trait_1   trait_2   trait_3
S1    S1    1.2       3.4       5.6
S2    S2    2.3       4.5       6.7
```

## Trait metadata format

Optional file with columns:

- `trait_id`
- `trait_level`
- `display_name`

If this file is not provided, the workflow automatically uses:

- `trait_id = column name`
- `trait_level = trait`
- `display_name = trait_id`

## Covariate format

The covariate file must contain:

- column 1: `FID`
- column 2: `IID`
- column 3 onward: one or more numeric covariates

## Environment variables

- `GWAS_WORKSPACE_DIR`
- `GWAS_DATA_DIR`
- `GWAS_TRAIT_MATRIX`
- `GWAS_TRAIT_META` (optional)
- `GWAS_TPED_PREFIX`
- `GWAS_KINSHIP`
- `GWAS_COVARIATES`
- `GWAS_FASTA_FAI_FILE` (optional)
- `GWAS_THREADS`
- `GWAS_EMMAX_PARALLEL`
- `EMMAX_BIN`

## Notes

- If `GWAS_FASTA_FAI_FILE` is missing, chromosome lengths are inferred from the genotype position file.
- This skill assumes `EMMAX`-ready inputs already exist.
