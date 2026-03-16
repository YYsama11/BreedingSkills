# Skill Index

## Standard Layout

- `SKILL.md`: trigger conditions, workflow, input requirements, and output expectations
- `scripts/`: reusable scripts for deterministic setup or execution
- `references/`: detailed materials that should be loaded only when needed
- `assets/`: templates and reusable resource files

## Current Skills

- `read-qc-alignment`: from raw FASTQ to sorted BAM and mapping QC
- `variant-call-filter`: from aligned BAM to normalized and filtered cohort VCF
- `genotype-kinship-prep`: from genotype plus phenotype to `emmax-intel64` inputs
- `gwas-run-visualization`: batch EMMAX execution plus Manhattan, QQ, and lead SNP outputs
- `qtl-candidate-gene`: convert significant SNPs into QTL regions and candidate gene tables
