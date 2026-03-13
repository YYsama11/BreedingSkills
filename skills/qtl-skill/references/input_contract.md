# Input Contract

## Required inputs

- `master_significant_snps.tsv.gz`
- PLINK genotype prefix:
  - `.bed`
  - `.bim`
  - `.fam`

## Optional inputs

- `annotation.gff3`
- `gene_annotation.tsv`
- `trait_summaries.tsv`

## Significant SNP table expectations

The GWAS significant SNP table should contain at least:

- `snp_id`
- `chrom`
- `pos`
- `p`
- `trait_id`

## Optional gene annotation table

If provided, this table should contain:

- `gene_id`

Optional extra columns may include:

- `description`
- `go`
- `kegg`
- `pfam`
- `interpro`

The workflow keeps all additional columns and merges them into candidate-gene outputs when possible.
