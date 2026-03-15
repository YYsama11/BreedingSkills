# Workflow

## Phase 1 — Lead locus selection

- `scripts/prepare_qtl_leads.py`

Input:

- GWAS significant SNP table

Output:

- `global_lead_snps.tsv`
- `global_lead_snps.txt`

## Phase 2 — LD calculation

- `scripts/run_ld_batch.py`
- `scripts/run_ld_for_lead.sh`

Input:

- global lead SNP list
- PLINK genotype prefix

Output:

- one LD file per lead locus

## Phase 3 — Annotation preparation (optional)

- `scripts/prepare_gene_annotation.py`

Input:

- GFF/GTF

Output:

- gene / promoter / exon / CDS interval resources

## Phase 4 — QTL integration

- `scripts/integrate_qtl_candidates.py`

Outputs:

- QTL intervals
- hotspot table
- trait membership table
- candidate genes
- representative genes

## Phase 5 — QTL plotting

- `scripts/build_qtl_plot_tables.py`
- `scripts/plot_gwas_qtl_summary.py`
- `scripts/plot_ld_qtl_summary.R`
- `scripts/run_qtl_plot.sh`

Inputs:

- either plotting tables described in `references/plot_input_contract.md`
- or QTL/GWAS outputs that can be converted into those plotting tables

Outputs:

- global Manhattan + local QTL + gene track summary figure

The plotting layer is intentionally decoupled from the core QTL logic so users can provide:

- only the minimal plotting tables
- or richer tables with highlight genes and LD values

Annotation-dependent outputs are skipped if annotation is unavailable.
