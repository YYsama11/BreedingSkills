# Workflow Notes

## 批量运行建议

- 先统一输出 trait 列表，再批量生成 phenotype 向量
- 每个性状一个结果前缀，避免覆盖
- 把日志与关联结果分目录存放，便于回查失败 trait

## 可视化最小集合

- Manhattan plot
- QQ plot
- lead SNP 表
- 每个 trait 的显著位点统计

## 交付标准

- 每个 trait 都能映射到唯一结果文件
- 结果表包含染色体、坐标、p 值和 SNP 标识
- 可视化输入列在所有 trait 间保持一致
