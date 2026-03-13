---
name: gwas-skill
description: Run and organize a high-throughput EMMAX-based GWAS workflow for many phenotypes, including replicate averaging, INT transformation, PLINK genotype preparation, EMMAX association, Manhattan/QQ summaries, global lead-locus QTL calling, candidate gene extraction from GFF, omics/network-style phenotype analyses, and subgroup allele-frequency comparisons. Use when a user wants to package, rerun, adapt, or publish a full large-scale GWAS workflow.
---

# GWAS Skill

Use this skill when the user wants a reusable GWAS project rather than one-off ad hoc commands.

## Quick start

1. Read `references/input_contract.md` to confirm required files, default filenames, and environment variables.
2. Read `references/workflow.md` for the end-to-end phase order and output layout.
3. For large runs, read `references/performance.md` before launching anything heavy.
4. Run the full workflow with:

```bash
bash scripts/run_full_pipeline.sh
```

## What this skill packages

- Replicate-aware phenotype generation with mean aggregation by sample prefix
- INT-transformed phenotype matrices for molecule/class/superclass/role/bin/total levels
- PLINK subsetting, PCA, allele frequencies, and kinship submatrix alignment
- EMMAX batch GWAS with per-trait logging and resumable downstream summarization
- Manhattan/QQ plotting and master significant SNP summaries
- Global lead-locus QTL calling using LD around nonredundant significant SNPs
- Candidate gene extraction from GFF-derived gene/promoter/CDS/exon tracks
- Lipid/network/population post-analysis built from GWAS outputs

## When to read the references

- **Input naming or directory changes**: `references/input_contract.md`
- **Pipeline order or output mapping**: `references/workflow.md`
- **Large data / speed / thread control**: `references/performance.md`
- **What happened in the reference GWAS run**: `references/observed_run_summary.md`

## Key operational rules

- Prefer environment variables over patching paths directly.
- Keep `GWAS_THREADS` and `GWAS_EMMAX_PARALLEL` conservative on shared servers.
- For many significant loci, use the global lead-locus workflow rather than per-trait LD expansion.
- Treat SNP functional consequence calls as unsupported unless real nucleotide alleles are available.
