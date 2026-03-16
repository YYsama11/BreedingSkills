# Input Contract

## 必需文件

- `PREFIX.tped`
- `PREFIX.tfam`
- kinship 矩阵
- `phenotype.tsv`

## 可选文件

- `covariates.tsv`

## `phenotype.tsv` 列要求

- 第 1 列必须为 `sample_id`
- 第 2 列开始为各性状名
- 性状名建议只使用字母、数字、下划线

## 结果表建议列

- `trait`
- `chrom`
- `pos`
- `snp_id`
- `pvalue`
- `beta`
- `stderr`
