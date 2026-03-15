# GWAS Skill

`gwas-skill` provides the **core GWAS analysis layer** of the repository.

It is intentionally limited to the stage that starts from:

- phenotype matrix
- `EMMAX`-ready genotype inputs

and ends with:

- Manhattan plots
- QQ plots
- significant SNP tables
- GWAS summary tables

This skill no longer bundles project-specific phenotype engineering, QTL interpretation, or population-specific downstream analysis into the default workflow.

---

## Scope

This skill is designed for users who already have:

- a phenotype matrix with `FID` and `IID`
- a genotype prefix in `tped/tfam` format
- a kinship matrix
- a covariate file

Optional input:

- chromosome lengths (`.fai`) for cleaner genome-wide plotting

If `FAI` is not provided, chromosome lengths are inferred from the maximum SNP position per chromosome in `TPED` or `BIM`.

---

## Required outputs

- Manhattan plot(s)
- QQ plot(s)

Additional core outputs:

- significant SNP table
- shared SNP table
- pleiotropic SNP table
- per-trait summary table

---

## Repository layout

```text
skills/gwas-skill/
├── SKILL.md
├── README.md
├── references/
│   ├── input_contract.md
│   ├── workflow.md
│   └── performance.md
└── scripts/
    ├── run_full_pipeline.sh
    ├── run_with_top.sh
    ├── prepare_emmax_inputs.py
    ├── run_emmax_batch.py
    ├── run_emmax_trait.sh
    ├── prepare_snp_reference.py
    ├── run_gwas_plot.sh
    ├── plot_manhattan_qq.R
    └── summarize_emmax_results.py
```

---

## Input model

### Mandatory inputs

- phenotype matrix
- `tped/tfam`
- kinship matrix
- covariate file

### Optional inputs

- trait metadata table
- chromosome length file (`FAI`)

Full details:

- `references/input_contract.md`

---

## Quick start

### 1. Prepare an external workspace

Example layout:

```text
your_project/
├── data/
│   ├── trait_matrix.tsv
│   ├── trait_metadata.tsv              # optional
│   ├── genotype_panel.tped
│   ├── genotype_panel.tfam
│   ├── kinship_matrix.tsv
│   ├── covariates.tsv
│   └── reference.fa.fai               # optional
└── analysis/
```

### 2. Set environment variables

```bash
export GWAS_WORKSPACE_DIR=/path/to/your_project
export GWAS_DATA_DIR=/path/to/your_project/data
export GWAS_TRAIT_MATRIX=/path/to/your_project/data/trait_matrix.tsv
export GWAS_TRAIT_META=/path/to/your_project/data/trait_metadata.tsv   # optional
export GWAS_TPED_PREFIX=/path/to/your_project/data/genotype_panel
export GWAS_KINSHIP=/path/to/your_project/data/kinship_matrix.tsv
export GWAS_COVARIATES=/path/to/your_project/data/covariates.tsv
export GWAS_FASTA_FAI_FILE=/path/to/your_project/data/reference.fa.fai # optional
export GWAS_THREADS=6
export GWAS_EMMAX_PARALLEL=10
```

### 3. Run the workflow

Run this command from inside `skills/gwas-skill/`:

```bash
bash scripts/run_full_pipeline.sh
```

`GWAS_WORKSPACE_DIR` is the external analysis workspace.  
The scripts remain inside the skill directory and are not expected to be copied into the workspace.

---

## Core outputs

### EMMAX raw results

- `analysis/gwas/emmax/results/*.ps`
- `analysis/gwas/emmax/results/*.reml`

### Summaries

- `analysis/gwas/emmax/trait_summaries.tsv`
- `analysis/gwas/emmax/master_significant_snps.tsv.gz`
- `analysis/gwas/emmax/shared_snp_stats.tsv`
- `analysis/gwas/emmax/pleiotropic_snps.tsv`
- `analysis/gwas/emmax/chromosome_distribution.tsv`

### Figures

- `analysis/gwas/emmax/figures/manhattan/*.png`
- `analysis/gwas/emmax/figures/qq/*.png`

Optional standalone plotting entrypoint:

- `scripts/run_gwas_plot.sh`

---

## Design constraints

- This skill assumes GWAS-ready genotype inputs and does not perform variant calling.
- This skill does not assume any phenotype domain such as lipids or metabolites.
- This skill does not require annotation or QTL support.
- QTL, LD, and candidate-gene interpretation are intentionally separated into `qtl-skill`.
- A standalone R plotting script is included for Manhattan + QQ rendering from generic GWAS tables.
- By default, tool-generated `.log` files are removed after successful execution.
- By default, runtime tracing keeps `analysis/logs/runtime/command_manifest.tsv` and does not retain per-command `*.cmd.log` / `*.top.log` files unless requested.

---

## Recommended companion skills

- `skills/resequencing-prep-skill` for generating EMMAX-ready inputs from raw reads
- `skills/qtl-skill` for LD, QTL, and candidate-gene interpretation
