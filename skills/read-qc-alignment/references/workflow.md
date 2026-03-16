# Workflow Notes

## 默认工具

- 原始质控：`fastqc`
- 过滤：`fastp`
- 比对：`bwa mem`
- 排序与索引：`samtools`

## 关键决策

- 如果是短读长群体重测序，优先使用 `bwa mem`
- 如果有固定接头序列或 lane 差异，先在样本表外部统一命名
- 如果后续要做 joint calling，尽量在本阶段统一 read group 命名规则

## 交付标准

- 每个样本至少有一个排序且已索引的 BAM
- 每个样本有 `flagstat` 或同等级统计文件
- 过滤前后均保留 QC 结果，便于回溯
