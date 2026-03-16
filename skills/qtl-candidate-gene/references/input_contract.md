# Input Contract

## 必需文件

- `gwas_summary.tsv`
- `gene_annotation.tsv`

## `gwas_summary.tsv` 列要求

- 表头至少包含：`trait	chrom	pos	snp_id	pvalue`

## `gene_annotation.tsv` 列要求

- 表头至少包含：`chrom	start	end	gene_id`

## 可选文件

- `ld_table.tsv`

## 输出建议

- QTL 区间表
- lead SNP 表
- candidate gene 表
