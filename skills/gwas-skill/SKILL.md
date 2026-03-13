---
name: gwas-skill
description: Run a core high-throughput EMMAX-based GWAS workflow starting from phenotype matrices plus EMMAX-ready genotype inputs (tped/tfam, kinship, covariates). Use when a user wants Manhattan plots, QQ plots, significant SNP tables, and GWAS summary tables without forcing QTL, annotation, or project-specific phenotype logic.
---

# GWAS Skill

Use this skill for the **core GWAS stage** only.

## Scope

This skill starts from:

- phenotype matrix
- `tped/tfam`
- kinship matrix
- covariate file

It produces:

- Manhattan plots
- QQ plots
- significant SNP summaries
- per-trait GWAS summary tables

## What this skill does not require

This skill does **not** require:

- raw resequencing reads
- genome annotation
- functional annotation tables
- LD / QTL analysis

Those belong in separate skills.

## Read before use

- `references/input_contract.md`
- `references/workflow.md`
- `references/performance.md`

## Main entrypoint

```bash
bash scripts/run_full_pipeline.sh
```
