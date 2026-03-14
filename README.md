# BreedingSkills

`BreedingSkills` is a repository of modular analysis skills designed for breeding workflows.

Each skill targets a specific stage of a breeding-oriented genomics pipeline, so users can start from the stage that matches their current data.

## Available skills

<details>
<summary><strong>resequencing-prep-skill</strong></summary>

### What it does

- starts from raw resequencing reads
- aligns reads to a reference genome
- calls SNPs
- prepares GWAS-ready genotype resources

### Input requirements

- sample manifest
- FASTQ files
- reference genome

### Output requirements

- `genotype_panel.tped`
- `genotype_panel.tfam`
- PLINK files
- kinship matrix
- PCA
- covariates

</details>

<details>
<summary><strong>gwas-skill</strong></summary>

### What it does

- runs the core GWAS stage with `EMMAX`
- summarizes association results
- produces Manhattan and QQ plots

### Input requirements

- phenotype matrix
- `EMMAX`-ready genotype inputs
- kinship matrix
- covariate file

### Output requirements

- Manhattan plots
- QQ plots
- significant SNP table
- GWAS summary table

</details>

<details>
<summary><strong>qtl-skill</strong></summary>

### What it does

- interprets GWAS signals at the locus level
- computes LD-supported QTL intervals
- annotates SNP positions
- extracts candidate genes
- renders QTL summary plots

### Input requirements

- significant SNP table from GWAS
- genotype files for LD calculation
- optional annotation files

### Output requirements

- lead loci
- QTL intervals
- hotspot table
- candidate gene table
- representative gene table
- QTL summary plots

</details>
