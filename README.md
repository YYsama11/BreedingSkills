# BreedingSkills

`BreedingSkills` is a modular collection of reusable breeding-analysis skills.

The repository is organized around separate analysis stages so that users can enter the workflow at the stage that matches their data:

- raw resequencing reads
- GWAS-ready genotype resources
- GWAS significant SNP results

---

## Repository structure

```text
BreedingSkills/
├── LICENSE
├── README.md
└── skills/
    ├── resequencing-prep-skill/
    ├── gwas-skill/
    ├── qtl-skill/
    ├── transposon-skill/     # future
    ├── pangenome-skill/      # future
    └── ...
```

---

## Included skills

### 1. `resequencing-prep-skill`

Purpose:

- start from raw resequencing reads
- align reads to a reference genome
- call SNPs
- convert results into `EMMAX`-ready genotype resources

Outputs include:

- `genotype_panel.tped/tfam`
- PLINK files
- kinship matrix
- PCA
- covariates

### 2. `gwas-skill`

Purpose:

- run the core GWAS stage
- start from phenotype matrix plus `EMMAX`-ready genotype inputs
- generate Manhattan plots, QQ plots, and significant SNP summaries

This skill intentionally excludes:

- raw read processing
- QTL interpretation
- project-specific phenotype subclassification

### 3. `qtl-skill`

Purpose:

- interpret GWAS signals at the locus level
- generate lead loci
- compute LD-supported QTL intervals
- annotate SNP positions
- collect candidate genes
- identify representative genes per QTL

---

## Design philosophy

This repository separates three different analysis responsibilities:

1. **resequencing preparation**
2. **core GWAS**
3. **QTL and candidate-gene interpretation**

This separation improves:

- portability
- maintainability
- reuse across projects
- compatibility with different data starting points

---

## Recommended use

- Start with `resequencing-prep-skill` if you only have raw `FASTQ` files.
- Start with `gwas-skill` if you already have `EMMAX`-ready inputs.
- Start with `qtl-skill` if you already have GWAS significant SNP outputs and want downstream biological interpretation.

---

## Planned future skills

- `transposon-skill`
- `pangenome-skill`
- `sv-skill`
- additional breeding-oriented genomics workflows
