# Skill Index

This repository is organized as a skill collection.

## Available skills

### `resequencing-prep-skill`

Use when the starting point is:

- raw resequencing reads
- reference genome
- sample manifest

Main outputs:

- GWAS-ready genotype resources
- kinship
- PCA
- covariates

### `gwas-skill`

Use when the starting point is:

- phenotype matrix
- `EMMAX`-ready genotype inputs

Main outputs:

- Manhattan plots
- QQ plots
- significant SNP summaries

### `qtl-skill`

Use when the starting point is:

- GWAS significant SNP results
- genotype files for LD calculation

Main outputs:

- QTL intervals
- candidate genes
- representative genes
- QTL summary plots

## Recommended order

Typical full workflow:

1. `resequencing-prep-skill`
2. `gwas-skill`
3. `qtl-skill`

Users do not need to run all three.  
Each skill is designed to be entered independently if the required inputs already exist.
