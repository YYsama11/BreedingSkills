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

## Use with Codex

This repository is designed to be:

- published on GitHub
- cloned locally
- installed into a local Codex skill directory

### Recommended usage model

1. Keep GitHub as the publication and update source
2. Clone the repository locally
3. Install one or more skills into your local Codex skills directory

### Install all skills

```bash
git clone https://github.com/YYsama11/BreedingSkills.git
cd BreedingSkills
bash scripts/install_to_codex_home.sh
```

### Install one skill only

```bash
bash scripts/install_to_codex_home.sh gwas-skill
```

### Install as symlinks for development

```bash
bash scripts/install_to_codex_home.sh --symlink
```

By default, the install target is:

- `${CODEX_HOME:-$HOME/.codex}/skills`

### Skill order

See:

- `skills/INDEX.md`

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
