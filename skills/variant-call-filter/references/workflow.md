# Workflow Notes

## caller 选择

- 样本数中等且追求快速迭代时，可先使用 `bcftools`
- 如果项目已经标准化到 GATK Best Practices，可切换到 `GATK HaplotypeCaller` + joint genotyping

## 过滤建议

- 先做位点标准化，再应用质量阈值
- 过滤至少包含 `QUAL`、深度、缺失率
- 群体数据常常还要追加 MAF 或 MAC 限制

## 交付标准

- 过滤后的 VCF 可以被 `plink` 或 `bcftools stats` 正常读取
- 位点数变化在各步骤中可解释
- 样本顺序与上游 BAM 清单可追溯
