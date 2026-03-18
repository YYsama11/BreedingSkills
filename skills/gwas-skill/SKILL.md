---
name: gwas-skill
description: Use when EMMAX-ready genotype, kinship, phenotype, and optional covariate files are already prepared and you need batch GWAS execution, result collection, and visualization inputs such as Manhattan, QQ, and lead SNP tables, or when you need help deciding whether to run GWAS on raw phenotypes, INT-transformed values, BLUEs, or BLUPs.
---

# gwas-skill

## When to use

- `tped/tfam`, kinship, and phenotype inputs are already prepared
- You need to run `emmax-intel64` across many traits
- You want a consistent set of Manhattan, QQ, and lead SNP outputs
- You need help deciding whether phenotype preprocessing should use raw values, INT, BLUE, or BLUP before GWAS

## Quick start

- Prepare EMMAX inputs and phenotype matrices described in `references/input_contract.md`
- Copy and edit `assets/gwas_run_visualization.env.template`
- Use `emmax-intel64` as the default association binary; do not probe plain `emmax` first unless your environment explicitly uses that name
- Before running GWAS, confirm with the user whether the phenotype file contains raw values, line means, BLUEs, BLUPs, or INT-transformed values
- Run `bash scripts/run_gwas_visualization.sh --config assets/gwas_run_visualization.env.template`

## Before running

- Do not silently decide between raw phenotype values, INT, BLUE, or BLUP
- Ask the user what phenotype scale is currently provided if it is not explicit
- Ask whether the data include replicates, blocks, years, or multiple environments
- Ask whether a mixed model has already been fitted upstream
- If the answer is unclear, pause execution and use `references/phenotype_decision.md` to guide the recommendation before launching EMMAX

## Default workflow

1. Confirm with the user whether GWAS should use raw values, INT, BLUEs, or BLUPs
2. Validate `tped/tfam`, kinship, phenotype, and optional covariate inputs
3. Build a trait list from the phenotype matrix
4. Generate per-trait phenotype vectors and EMMAX commands
5. Collect association results, significant hits, and plotting tables
6. Export a standard directory structure for Manhattan, QQ, and summary figures

## Primary outputs

- `trait_list.txt`
- `traits/`
- `assoc/`
- `run_emmax_commands.sh`
- `plot_plan.txt`

## Read only when needed

- Input file and column requirements: `references/input_contract.md`
- Result directory organization guidance: `references/workflow.md`
- Plot colors, panel layout, and threshold style: `references/plot_style.md`
- Phenotype decision guide for INT, BLUE, and BLUP: `references/phenotype_decision.md`
