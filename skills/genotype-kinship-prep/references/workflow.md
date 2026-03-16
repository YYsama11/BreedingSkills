# Workflow Notes

## 默认软件

- 格式转换与 PCA：`plink`
- kinship：`emmax-kin-intel64`
- 关联分析引擎：`emmax-intel64`

## 关键检查点

- phenotype 首列样本顺序是否与 `tfam` 对齐
- PCA 与 kinship 是否基于同一批样本
- covariate 列是否包含截距或已按项目规范处理

## 交付标准

- 生成的 `tped/tfam` 可直接供 `emmax-intel64` 调用
- kinship 文件与 `tfam` 样本顺序一致
- 所有 trait 已汇总到一个清晰的 trait 列表中
