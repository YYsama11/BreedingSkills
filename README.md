# BreedingSkill

`BreedingSkill` is a growing repository of reusable breeding-analysis skills.

This project is designed as a **multi-skill collection**, where each skill lives in its own subdirectory under `skills/` and can be maintained independently.

Current and future directions include:

- GWAS
- transposon analysis
- pangenome analysis
- structural variation analysis
- population genetics
- other breeding-oriented genomics workflows

## Repository layout

```text
BreedingSkill/
├── LICENSE
├── README.md
└── skills/
    ├── gwas-skill/
    ├── transposon-skill/     # future
    ├── pangenome-skill/      # future
    └── ...
```

## Included skills

### `gwas-skill`

Location:

- `skills/gwas-skill`

This skill packages a reusable high-throughput GWAS workflow derived from a real large-scale omics association run. It includes:

- replicate-aware phenotype aggregation
- INT transformation
- PLINK-based genotype preparation
- EMMAX batch GWAS
- Manhattan / QQ / significant SNP summaries
- global lead-locus QTL calling
- candidate gene extraction
- network and subgroup frequency analysis

See:

- `skills/gwas-skill/README.md`
- `skills/gwas-skill/SKILL.md`

## Expansion strategy

New skills should be added as siblings under `skills/`, for example:

- `skills/transposon-skill`
- `skills/pangenome-skill`
- `skills/sv-skill`

This keeps the repository modular and suitable both as a GitHub project and as a reusable analysis skill collection.
