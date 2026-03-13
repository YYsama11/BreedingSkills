# QTL Skill

`qtl-skill` is the downstream interpretation layer that turns GWAS signals into locus-level outputs.

It is designed to start **after GWAS is complete**.

---

## Scope

Mandatory inputs:

- significant SNP table from GWAS
- PLINK-compatible genotype files for LD calculation

Optional inputs:

- annotation GFF/GTF
- gene functional annotation table
- trait summary table

Mandatory outputs:

- lead loci
- LD-based QTL intervals
- QTL hotspot table

Optional outputs:

- SNP positional annotation
- all genes inside each QTL
- representative gene per QTL

SNP positional annotation categories include:

- intergenic
- upstream 2 kb
- intron
- exon
- downstream 2 kb

---

## Input model

### Mandatory

- `master_significant_snps.tsv.gz` or equivalent
- PLINK `BED/BIM/FAM` prefix for LD calculation

### Optional

- `annotation.gff3`
- `gene_annotation.tsv`
- `trait_summaries.tsv`

If annotation is missing, the skill still produces:

- lead loci
- LD intervals
- QTL regions
- hotspot summaries

but skips gene-based outputs.

---

## Main outputs

- `analysis/qtl/global_lead_snps.tsv`
- `analysis/qtl/qtl_regions.tsv`
- `analysis/qtl/qtl_hotspots.tsv`
- `analysis/qtl/qtl_trait_membership.tsv`
- `analysis/qtl/qtl_trait_counts.tsv`
- `analysis/qtl/qtl_size_stats.tsv`

If annotation is available:

- `analysis/qtl/significant_snp_annotation.tsv.gz`
- `analysis/qtl/candidate_genes.tsv`
- `analysis/qtl/qtl_representative_genes.tsv`

---

## Main scripts

```text
scripts/
├── run_qtl_pipeline.sh
├── run_with_top.sh
├── prepare_gene_annotation.py
├── prepare_qtl_leads.py
├── run_ld_batch.py
├── run_ld_for_lead.sh
└── integrate_qtl_candidates.py
```

---

## Notes

- Hotspots are built from merged overlapping QTL intervals, not from a simple filter on large trait counts.
- Representative genes are workflow heuristics, not confirmed causal genes.
- If annotation chromosome names and GWAS chromosome names differ, the workflow attempts normalized matching rather than forcing `ChrN` naming.
