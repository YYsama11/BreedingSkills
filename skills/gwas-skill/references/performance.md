# Performance Notes

## Main scaling factors

Core GWAS runtime is driven mainly by:

- number of traits
- number of SNPs
- degree of trait-level parallelism

## Current design choices

### 1. Trait-wise EMMAX execution

`EMMAX` is used trait by trait in parallel.  
This is simple, explicit, and works well for large phenotype panels.

### 2. Resumable summary generation

The GWAS summarization script caches:

- per-trait summaries
- top-hit tables
- significant SNP tables
- Manhattan plots
- QQ plots

This avoids recomputing all outputs after interruption.

### 3. Conservative thread wrapper

`scripts/run_with_top.sh`:

- checks `top`
- estimates a safe thread count
- exports `GWAS_THREADS`
- logs every command into `analysis/logs/runtime/command_manifest.tsv`

## Practical advice

- Keep `GWAS_EMMAX_PARALLEL` modest on shared servers
- Keep `GWAS_THREADS` modest when combining multiple parallel layers
- Use the optional `FAI` input if you want cleaner cumulative chromosome plots
