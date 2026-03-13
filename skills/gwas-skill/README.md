# GWAS Skill

`gwas-skill` is a reusable high-throughput GWAS workflow package derived from a fully executed large-scale omics association study.

The goal of this project is not to preserve a single one-off analysis, but to provide a workflow that is:

- reusable
- parameterized
- resumable
- suitable for large phenotype panels
- structured as a maintainable analysis skill

---

## Scope

This skill is designed for projects that need:

- batch GWAS over many phenotypes
- technical replicate aggregation before GWAS
- inverse normal transformed phenotypes
- mixed-model GWAS with `EMMAX`
- downstream QTL, candidate-gene, network, and subgroup frequency summaries

Typical use cases include metabolite GWAS, lipid GWAS, ionomics GWAS, and other high-dimensional quantitative phenotype studies.

---

## What this workflow solves

Large GWAS projects often fail at the workflow level rather than at the single-command level. Common problems include:

- too many input files to keep aligned correctly
- too many phenotypes to run manually
- output files becoming too large for practical post-processing
- QTL expansion becoming computationally unrealistic
- server-side overuse of threads and I/O

This repository addresses those issues by:

- standardizing inputs with environment variables
- tracking runtime decisions with `run_with_top.sh`
- using `EMMAX` for scalable trait-wise association scans
- caching downstream summaries so reruns are cheap
- collapsing significant loci into global nonredundant lead loci before LD/QTL steps
- reducing frequency summaries to SNP-level outputs instead of trait-SNP-level explosion

---

## Workflow overview

The full workflow is organized into these phases:

1. **Phenotype preprocessing**
   - detect three replicate columns per sample
   - average replicates by sample prefix
   - build molecule, class, superclass, role, chain-bin, and total-abundance traits
   - apply inverse normal transformation

2. **Annotation resource preparation**
   - derive gene, promoter, exon, and CDS interval tables from GFF
   - backfill gene-level functional annotations from transcript features when needed

3. **Genotype preprocessing**
   - subset the genotype panel to target accessions
   - reorder covariates and kinship consistently
   - build BED, TPED, PCA, and frequency resources

4. **EMMAX input preparation**
   - write one phenotype file per trait
   - create a trait manifest

5. **Batch GWAS with EMMAX**
   - run trait-wise mixed-model association in parallel

6. **GWAS summarization**
   - Manhattan plots
   - QQ plots
   - significant SNP tables
   - shared SNP and pleiotropy summaries

7. **QTL and LD analysis**
   - collapse significant SNPs into global lead loci
   - compute LD around each lead locus
   - define QTL intervals and hotspots

8. **Candidate gene analysis**
   - extract genes overlapping QTL intervals
   - mark genes with keyword-based functional categories

9. **Network and population analysis**
   - phenotype structure analysis
   - lipid or metabolite correlation networks
   - module assignments
   - subgroup allele-frequency summaries
   - selection-proxy overlap with QTL

---

## Repository layout

```text
gwas-skill/
├── SKILL.md
├── README.md
├── .gitignore
├── references/
│   ├── input_contract.md
│   ├── workflow.md
│   ├── performance.md
│   └── observed_run_summary.md
└── scripts/
    ├── run_full_pipeline.sh
    ├── run_post_emmax_pipeline.sh
    ├── run_with_top.sh
    ├── prepare_phenotypes.py
    ├── prepare_gene_annotation.py
    ├── prepare_genotypes.sh
    ├── prepare_emmax_inputs.py
    ├── run_emmax_batch.py
    ├── run_emmax_trait.sh
    ├── summarize_emmax_results.py
    ├── prepare_qtl_leads.py
    ├── run_ld_batch.py
    ├── run_ld_for_lead.sh
    ├── integrate_qtl_candidates.py
    ├── phenotype_network_population_analysis.py
    └── summarize_population_frequencies.py
```

---

## Quick start

### 1. Prepare the data directory

The default layout expects:

```text
your_project/
├── data/
│   ├── sample_ids.txt
│   ├── phenotype_matrix.tsv
│   ├── genotype_panel.tped
│   ├── genotype_panel.tfam
│   ├── kinship_matrix.tsv
│   ├── covariates.tsv
│   ├── annotation.gff3
│   └── reference.fa.fai
└── analysis/
```

Full input expectations are documented in:

- `references/input_contract.md`

### 2. Set environment variables

Minimal example:

```bash
export GWAS_WORKSPACE_DIR=/path/to/your_project
export GWAS_DATA_DIR=/path/to/your_project/data
export GWAS_THREADS=6
export GWAS_EMMAX_PARALLEL=10
```

`GWAS_WORKSPACE_DIR` is the external analysis workspace that contains `data/` and will receive `analysis/`.
The workflow scripts themselves are executed from the `skills/gwas-skill/` repository directory and do not need to be copied into the workspace.

If tools are not in `PATH`, define them explicitly:

```bash
export PLINK_BIN=/path/to/plink
export EMMAX_BIN=/path/to/emmax-intel64
```

### 3. Run the full workflow

```bash
bash scripts/run_full_pipeline.sh
```

---

## Main inputs

### Phenotypes

- `sample_ids.txt`: target sample identifiers
- `phenotype_matrix.tsv`: phenotype matrix with metadata columns followed by replicate columns

Assumptions:

- replicate columns use names such as `Sample_1`, `Sample_2`, `Sample_3`
- replicate columns do not need to be adjacent
- grouping is based on the sample prefix, not on fixed column position

### Genotypes

- `*.tped` / `*.tfam`
- kinship matrix
- covariate table

### Annotation

- `GFF`
- `FAI`

---

## Main outputs

### GWAS summaries

- `analysis/gwas/emmax/trait_summaries.tsv`
- `analysis/gwas/emmax/master_significant_snps.tsv.gz`
- `analysis/gwas/emmax/shared_snp_stats.tsv`
- `analysis/gwas/emmax/figures/manhattan/`
- `analysis/gwas/emmax/figures/qq/`

### QTL and candidate genes

- `analysis/qtl/qtl_regions.tsv`
- `analysis/qtl/qtl_hotspots.tsv`
- `analysis/qtl/candidate_genes.tsv`
- `analysis/qtl/candidate_genes_lipid_related.tsv`

### Network and population outputs

- `analysis/network_population/lipid_network_nodes.tsv`
- `analysis/network_population/lipid_module_membership.tsv`
- `analysis/population_genetics/significant_snp_allele_frequencies.tsv`

---

## Performance strategy

This workflow has already been optimized against a real large-scale run. Key optimizations include:

### 1. EMMAX-based batch GWAS

- uses `EMMAX` for large numbers of univariate traits
- supports parallel trait execution

### 2. Resumable downstream summaries

- Manhattan, QQ, top-hit, and per-trait summary outputs are cached
- reruns do not regenerate completed outputs unnecessarily

### 3. Global lead-locus QTL design

Instead of expanding LD independently for every phenotype:

- significant SNPs are deduplicated
- global nonredundant lead loci are selected
- LD and QTL steps are run per lead locus

This keeps QTL processing tractable even when many traits share the same association peak.

### 4. Vectorized SNP annotation

- avoids slow row-wise interval scans
- uses chromosome-wise sorted interval lookups

### 5. SNP-level population summaries

- avoids trait-SNP-scale frequency tables
- writes SNP-level subgroup frequency outputs instead

More detail is documented in:

- `references/performance.md`

---

## Known limitations

### 1. Non-synonymous SNP calling requires real nucleotide alleles

If genotype alleles are encoded as `1/2` instead of `A/C/G/T`, codon-level consequence inference is not reliable.

### 2. Hotspot definitions are threshold-sensitive

If many phenotypes share the same large signal block, hotspot counts and sizes can change substantially with different thresholds.

### 3. Some QTL intervals may collapse to lead-only loci

Under the current LD thresholds, some loci may produce:

- `qtl_start == qtl_end`

This is expected behavior for sparse or weakly linked surrounding regions.

---

## Reference run scale

This workflow was extracted from a fully executed large-scale omics GWAS project with approximately:

- 112 accessions
- 3 replicates per accession
- 1461 molecular phenotypes
- 1534 total traits
- 4,033,947 SNPs

That run produced:

- 348 traits with significant signals
- 89,802 unique significant SNPs
- 552 global lead or QTL loci
- 621 unique candidate genes

See:

- `references/observed_run_summary.md`

---

## Recommended next additions

Useful future improvements for this skill include:

- a more general configuration file format
- stricter hotspot definitions for publication-grade reporting
- nucleotide-level consequence analysis when nucleotide alleles are available
- automated compression and archival of raw GWAS outputs
