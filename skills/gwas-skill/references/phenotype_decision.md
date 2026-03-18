# Phenotype Decision Guide

Use this guide before launching GWAS if the user has not already specified whether the phenotype should be analyzed as raw values, INT-transformed values, BLUEs, or BLUPs.

## Ask the user first

Ask these questions before running `emmax-intel64`:

- Are the current phenotype values raw plot-level observations, genotype-level means, BLUEs, BLUPs, or already INT-transformed values?
- Does the experiment include replicates, blocks, years, sites, or other structured trial effects?
- Has a mixed model already been fitted upstream?
- Does the user want GWAS on observed means, adjusted means, BLUPs, or transformed values?

Do not infer these choices silently when the answer matters to the model interpretation.

## Decision rules

### Use raw values directly when

- Each genotype has one final value already summarized at the genotype level
- There is no replicate/block/site/year structure left to model
- The phenotype distribution is acceptable for direct GWAS

### Consider INT when

- The phenotype is strongly skewed, heavy-tailed, or dominated by outliers
- The user explicitly wants rank-based normalization for GWAS
- The trait has already been reduced to one value per genotype

Avoid recommending INT as an automatic default for every trait.

### Consider BLUE or BLUP when

- The phenotype still comes from replicated or multi-environment trials
- Block, year, location, or other design effects should be modeled before GWAS
- The user wants adjusted genotype-level values rather than raw plot observations

### BLUE vs BLUP

- BLUE is more natural when the genotype effect is treated as fixed and the user wants adjusted means within the observed experiment
- BLUP is more natural when genotype effects are treated as random and the user wants shrinkage/prediction across a structured trial design

## Practical recommendation

- If the phenotype file is already one row per genotype and the user confirms it is final, run GWAS directly unless they request INT
- If the phenotype file still reflects trial structure, stop and recommend fitting a mixed model first to derive BLUEs or BLUPs
- If the user is unsure, recommend comparing at least two phenotype representations rather than pretending there is only one correct choice

## What this skill currently assumes

The current `gwas-skill` runner expects a finalized genotype-level phenotype matrix. It does not yet fit mixed models for BLUE/BLUP generation by itself. That choice must be made and prepared before launching the GWAS step.
