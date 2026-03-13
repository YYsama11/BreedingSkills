---
name: resequencing-prep-skill
description: Prepare GWAS-ready genotype inputs from raw resequencing reads by aligning FASTQ files to a reference genome, calling SNPs, converting to PLINK and EMMAX-compatible formats, building kinship, computing PCA, and generating covariates.
---

# Resequencing Preparation Skill

Use this skill when the user only has raw resequencing reads and wants to generate inputs required by the GWAS skill.

## Scope

This skill starts from:

- FASTQ files
- reference genome
- sample manifest

It produces:

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
