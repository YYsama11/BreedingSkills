---
name: qtl-skill
description: Interpret GWAS signals at the locus level by collapsing significant SNPs into lead loci, computing LD-supported QTL regions, annotating SNP positions relative to genes, collecting candidate genes, and selecting representative genes per QTL when annotation is available.
---

# QTL Skill

Use this skill after GWAS has already produced significant SNP results.

## Scope

This skill starts from:

- significant SNP tables
- genotype files for LD calculation

Optional inputs:

- GFF/GTF
- gene functional annotation table
- trait summary table

It produces:

- lead loci
- LD-supported QTL intervals
- hotspot summaries
- candidate gene tables
- representative gene tables

## Read before use

- `references/input_contract.md`
- `references/workflow.md`

## Main entrypoint

```bash
bash scripts/run_qtl_pipeline.sh
```
