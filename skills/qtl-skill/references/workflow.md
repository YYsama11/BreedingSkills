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

Annotation-dependent outputs are skipped if annotation is unavailable.
