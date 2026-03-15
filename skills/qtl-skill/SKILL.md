---
name: qtl-skill
description: Interpret GWAS signals at the locus level by collapsing significant SNPs into lead loci, computing LD-supported QTL regions, annotating SNP positions relative to genes, collecting candidate genes, selecting representative genes per QTL when annotation is available, and rendering layered LD/QTL summary figures from minimal or rich plotting tables.
---

# QTL Skill

Use this skill after GWAS has already produced significant SNP results.

## When to use

- The user already has GWAS significant SNP results.
- The user wants lead loci, LD-supported QTL regions, or candidate-gene interpretation.
- The user may also want layered QTL summary figures.

## When not to use

- Do not use this skill for raw read processing.
- Do not use this skill for core GWAS model execution.
- Use `gwas-skill` first if association analysis has not been run yet.

## Inputs

- significant SNP tables
- genotype files for LD calculation

Optional inputs:

- GFF/GTF
- gene functional annotation table
- trait summary table

## Outputs

- lead loci
- LD-supported QTL intervals
- hotspot summaries
- candidate gene tables
- representative gene tables
- layered LD/QTL summary figures

## Read before use

- `references/input_contract.md`
- `references/workflow.md`
- `references/plot_input_contract.md`

## Main entrypoints

```bash
bash scripts/run_qtl_pipeline.sh
```

For plotting-only usage:

```bash
bash scripts/run_qtl_plot.sh
```
