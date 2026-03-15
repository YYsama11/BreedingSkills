---
name: gwas-skill
description: Run a core high-throughput EMMAX-based GWAS workflow starting from phenotype matrices plus EMMAX-ready genotype inputs (tped/tfam, kinship, covariates). Use when a user wants Manhattan plots, QQ plots, significant SNP tables, and GWAS summary tables without forcing QTL, annotation, or project-specific phenotype logic.
---

# GWAS Skill

Use this skill for the core GWAS stage only.

## When to use

- The user already has `EMMAX`-ready genotype inputs.
- The user has a phenotype matrix or preformatted trait table.
- The user wants GWAS association results, Manhattan plots, and QQ plots.

## When not to use

- Do not use this skill if the user only has raw resequencing reads.
- Use `resequencing-prep-skill` first in that case.
- Do not use this skill for LD/QTL/candidate-gene interpretation.
- Use `qtl-skill` after GWAS if downstream locus interpretation is needed.

## Inputs

- phenotype matrix
- `tped/tfam`
- kinship matrix
- covariate file

## Outputs

- Manhattan plots
- QQ plots
- significant SNP summaries
- per-trait GWAS summary tables

Optional standalone plotting:

- Manhattan + QQ from a generic GWAS result table through `scripts/run_gwas_plot.sh`

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

## Downstream handoff

The outputs of this skill are intended to feed directly into:

- `qtl-skill`
