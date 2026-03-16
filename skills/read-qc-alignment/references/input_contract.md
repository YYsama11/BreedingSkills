# Input Contract

## 必需文件

- `sample_sheet.tsv`
- 参考基因组 `genome.fa`
- 原始双端 FASTQ 文件目录

## `sample_sheet.tsv` 列要求

- 表头必须为：`sample_id	read1	read2`
- `read1` 和 `read2` 填写相对 `FASTQ_DIR` 的文件名或绝对路径
- `sample_id` 不能重复

## 最低检查项

- 样本名在所有中间文件中保持一致
- FASTQ 文件对数与样本表行数一致
- 参考基因组已经准备好 `bwa` 与 `samtools faidx` 索引

## 推荐输出组织

- `outputs/read_qc_alignment/qc_raw`
- `outputs/read_qc_alignment/qc_trimmed`
- `outputs/read_qc_alignment/bam`
- `outputs/read_qc_alignment/stats`
