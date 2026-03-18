# Input Contract

## Required files

- `PREFIX.tped`
- `PREFIX.tfam`
- A kinship matrix
- `phenotype.tsv`

## Optional files

- `covariates.tsv`

## `phenotype.tsv` columns

- The first column must be `sample_id`
- Columns from the second onward are trait names
- Trait names should use letters, numbers, and underscores when possible
- The phenotype file should already represent the final genotype-level scale chosen for GWAS, such as raw genotype means, BLUEs, BLUPs, or INT-transformed values
- If that status is not explicit, confirm it with the user before running GWAS

## Recommended result columns

- `trait`
- `chrom`
- `pos`
- `snp_id`
- `pvalue`
- `beta`
- `stderr`
