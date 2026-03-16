# BreedingSkills

[English](README.md) | [简体中文](README.zh-CN.md)

BreedingSkills is a skill repository for breeding-oriented resequencing, GWAS, and QTL interpretation workflows.

The current version provides five core skills, and each skill follows the same layout:

- `SKILL.md`: required core instructions
- `scripts/`: optional example scripts
- `references/`: optional reference material
- `assets/`: optional templates

## Skills

- `skills/read-qc-alignment`: raw FASTQ quality control, trimming, and read alignment
- `skills/variant-call-filter`: variant calling, normalization, and filtering
- `skills/genotype-kinship-prep`: EMMAX input preparation, covariates, and kinship
- `skills/gwas-run-visualization`: EMMAX GWAS execution and result organization
- `skills/qtl-candidate-gene`: QTL interpretation and candidate gene collection

## Entry Points

- Skill index: `skills/INDEX.md`
- Per-skill instructions: each skill folder contains its own `SKILL.md`
