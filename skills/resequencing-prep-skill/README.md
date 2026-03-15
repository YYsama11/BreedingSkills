# Resequencing Preparation Skill

`resequencing-prep-skill` is the upstream genotype-preparation layer of the repository.

It is designed for users who start from raw resequencing reads and need to generate `EMMAX`-ready genotype resources.

---

## Scope

Inputs:

- sample manifest
- FASTQ files
- reference genome

Outputs:

- `genotype_panel.tped`
- `genotype_panel.tfam`
- `genotype_panel.bed`
- `genotype_panel.bim`
- `genotype_panel.fam`
- `kinship_matrix.tsv`
- `covariates.tsv`
- `sample_ids.txt`

By default:

- `sample_ids.txt` is kept as an interface file for downstream steps
- tool-generated `.log` files are removed after successful execution
- `.nosex` files are removed after successful execution

---

## Main stages

1. align reads to the reference genome
2. sort and index BAM files
3. jointly call SNPs
4. filter to biallelic SNPs
5. convert to PLINK and TPED/TFAM
6. compute kinship
7. compute PCA
8. derive covariates from PCA

---

## Main scripts

```text
scripts/
├── run_with_top.sh
├── run_resequencing_prep.sh
└── make_covariates_from_pca.py
```

---

## Notes

- This is an upstream preparation skill, not the GWAS engine itself.
- The resulting outputs are intended to feed directly into `skills/gwas-skill`.
- Tool choices are intentionally conservative and based on widely available utilities:
  - `bwa` or `bwa-mem2`
  - `samtools`
  - `bcftools`
  - `plink`
  - `emmax-kin-intel64`
