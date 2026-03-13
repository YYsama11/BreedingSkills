# Performance Notes

## Main bottlenecks

### 1. EMMAX over many traits

For thousands of phenotypes, the dominant cost is:

- number of phenotypes
- number of SNPs
- per-trait repeated mixed-model fitting

Recommended controls:

- `GWAS_EMMAX_PARALLEL=10` was workable in the reference run
- keep `GWAS_THREADS` modest on shared servers because PLINK and downstream array work can oversubscribe cores

### 2. GWAS summary generation

Naively re-reading and re-plotting every trait is slow. This skill keeps:

- per-trait cached summary files
- existing Manhattan/QQ/top-hit outputs
- resumable aggregation behavior

### 3. QTL calling

Do **not** expand LD per trait when many traits share the same locus.

Use:

- unique significant SNPs
- then global nonredundant lead loci
- then LD per global lead

This reduces impossible workloads to a tractable number of LD jobs.

### 4. SNP annotation

Avoid row-wise interval scans over all genes.

This skill uses:

- chromosome-wise sorting
- vectorized `searchsorted`
- unique significant SNPs instead of trait-SNP rows

## Practical thread policy

The bundled `scripts/run_with_top.sh`:

- snapshots `top` before each command
- chooses a conservative thread count
- exports `GWAS_THREADS`
- records commands in `analysis/logs/runtime/command_manifest.tsv`

## Recommended output discipline

- keep raw EMMAX result files until summaries finish
- only compress or clean `.ps` files after downstream summaries succeed
- avoid writing trait-SNP-scale tables when SNP-scale tables carry the same information
