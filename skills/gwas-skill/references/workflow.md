# Workflow

## Phase 1 — Phenotypes

- `scripts/prepare_phenotypes.py`
- Builds raw and INT phenotype matrices
- Generates molecule, class, superclass, role, chain-length, unsaturation, chain-count, and total-abundance traits

Outputs:

- `analysis/phenotypes/all_traits_raw.tsv`
- `analysis/phenotypes/all_traits_int.tsv`
- `analysis/phenotypes/trait_metadata.tsv`

## Phase 2 — Annotation

- `scripts/prepare_gene_annotation.py`
- Converts GFF into gene, promoter, exon, and CDS interval resources
- Backfills gene annotations from transcript records when needed

Outputs:

- `analysis/annotation/gene_metadata.tsv`
- `analysis/annotation/genes.bed`
- `analysis/annotation/promoters_2kb.bed`

## Phase 3 — Genotypes

- `scripts/prepare_genotypes.sh`
- Subsets the full genotype panel to target IDs
- Aligns covariates and kinship to subset order
- Creates PLINK BED and TPED representations
- Computes PCA and allele frequencies

Outputs:

- `analysis/genotype/target_subset.*`
- `analysis/genotype/target_subset_tped.*`
- `analysis/genotype/kinship_subset.tsv`
- `analysis/genotype/target_subset_pca.eigenvec`

## Phase 4 — EMMAX input preparation

- `scripts/prepare_emmax_inputs.py`
- Writes one phenotype file per trait plus a manifest

Outputs:

- `analysis/emmax_inputs/trait_manifest.tsv`
- `analysis/emmax_inputs/phenotypes/*.tsv`

## Phase 5 — EMMAX GWAS

- `scripts/run_emmax_batch.py`
- Parallel trait-wise EMMAX association

Outputs:

- `analysis/gwas/emmax/results/*.ps`
- `analysis/gwas/emmax/results/*.reml`
- `analysis/gwas/emmax/results/*.log`

## Phase 6 — GWAS summarization

- `scripts/summarize_emmax_results.py`
- Generates per-trait summary rows, top hits, significant hits, Manhattan plots, and QQ plots

Outputs:

- `analysis/gwas/emmax/trait_summaries.tsv`
- `analysis/gwas/emmax/master_significant_snps.tsv.gz`
- `analysis/gwas/emmax/shared_snp_stats.tsv`
- `analysis/gwas/emmax/figures/manhattan/*.png`
- `analysis/gwas/emmax/figures/qq/*.png`

## Phase 7 — Network and phenotype structure

- `scripts/phenotype_network_population_analysis.py`
- Builds subgroup assignments, aggregate PCA, clustering, correlation networks, and module membership

Outputs:

- `analysis/network_population/sample_groups.tsv`
- `analysis/network_population/lipid_network_nodes.tsv`
- `analysis/network_population/lipid_module_membership.tsv`

## Phase 8 — QTL and candidates

- `scripts/prepare_qtl_leads.py`
- `scripts/run_ld_batch.py`
- `scripts/integrate_qtl_candidates.py`

This workflow uses **global nonredundant lead SNPs** instead of per-trait lead SNPs to keep LD/QTL processing tractable for many traits.

Outputs:

- `analysis/qtl/global_lead_snps.tsv`
- `analysis/qtl/qtl_regions.tsv`
- `analysis/qtl/qtl_hotspots.tsv`
- `analysis/qtl/candidate_genes.tsv`

## Phase 9 — Population genetics summaries

- subgroup frequencies with PLINK
- `scripts/summarize_population_frequencies.py`

Outputs:

- `analysis/population_genetics/significant_snp_allele_frequencies.tsv`
- `analysis/population_genetics/selection_proxy_snps.tsv`
- `analysis/population_genetics/selection_proxy_qtl_overlap.tsv`
