# Input Contract

## 必需文件

- 群体 VCF 或已过滤的 PLINK 数据
- `phenotype.tsv`

## `phenotype.tsv` 列要求

- 第 1 列必须为 `sample_id`
- 第 2 列开始为一个或多个性状
- 缺失值建议统一为 `NA`

## 样本一致性要求

- genotype 与 phenotype 使用同一套样本命名
- 去除重复样本后再计算 PCA 和 kinship
- 如有群体结构协变量，建议单独提供一个 covariate 表

## EMMAX 输出要求

- `PREFIX.tped`
- `PREFIX.tfam`
- kinship 矩阵
- 每个性状可单独抽取的 phenotype 向量
