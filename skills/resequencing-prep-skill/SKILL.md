---
name: resequencing-prep-skill
description: Prepare GWAS-ready genotype inputs from raw resequencing reads by aligning FASTQ files to a reference genome, calling SNPs, converting to PLINK and EMMAX-compatible formats, building kinship, computing PCA, and generating covariates.
---

# Resequencing Preparation Skill

Use this skill when the user starts from raw resequencing reads and wants outputs that can be passed into the GWAS skill.

## When to use

- The user has `FASTQ` files, a reference genome, and a sample manifest.
- The user does not yet have `tped/tfam`, kinship, PCA, or covariates.
- The user wants to prepare `EMMAX`-ready genotype inputs before running GWAS.

## When not to use

- Do not use this skill if the user already has `tped/tfam`, kinship, and covariates.
- In that case, use `gwas-skill` directly.

## Inputs

- `FASTQ` files
- reference genome
- sample manifest

## Outputs

- genotype panel in `tped/tfam`
- PLINK `BED/BIM/FAM`
- kinship matrix
- PCA
- covariates
- sample ID list

## Read before use

- `references/input_contract.md`
- `references/workflow.md`

## Main entrypoint

```bash
bash scripts/run_resequencing_prep.sh
```

## Downstream handoff

The outputs of this skill are intended to feed directly into:

- `gwas-skill`
