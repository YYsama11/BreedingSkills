# Workflow

## Phase 1 — Reference indexing

- create `bwa` index if missing
- create `samtools faidx` index if missing

## Phase 2 — Alignment

- align each sample to the reference genome
- produce sorted BAM files
- produce BAM indices

## Phase 3 — Joint SNP calling

- build a BAM list
- run `bcftools mpileup`
- run `bcftools call`
- filter to biallelic SNPs

## Phase 4 — PLINK / EMMAX conversion

- generate `BED/BIM/FAM`
- generate `TPED/TFAM`
- compute kinship
- compute PCA
- derive covariates from top PCs

## Final outputs

- genotype panel
- kinship matrix
- PCA
- covariates
- sample ID list
