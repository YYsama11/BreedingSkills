---
name: gwas-run-visualization
description: Use when EMMAX-ready genotype, kinship, phenotype, and optional covariate files are already prepared and you need batch GWAS execution, result collection, and visualization inputs such as Manhattan, QQ, and lead SNP tables; 适用于使用 emmax-intel64 完成 GWAS 并整理可视化结果的场景。
---

# gwas-run-visualization

## 何时使用

- `tped/tfam`、kinship、phenotype 已经准备完成
- 需要对多个性状批量运行 `emmax-intel64`
- 需要统一整理 Manhattan、QQ 和 lead SNP 表

## 快速开始

- 按 `references/input_contract.md` 准备 EMMAX 输入和表型矩阵
- 复制并修改 `assets/gwas_run_visualization.env.template`
- 运行 `bash scripts/run_gwas_visualization.sh --config assets/gwas_run_visualization.env.template`

## 默认工作流

1. 校验 `tped/tfam`、kinship、phenotype 与 covariate 输入
2. 从 phenotype 矩阵生成 trait 列表
3. 为每个 trait 生成单独 phenotype 向量和 EMMAX 命令
4. 统一收集关联结果、显著位点和可视化输入表
5. 输出供 Manhattan / QQ / summary figure 使用的标准目录

## 主要产物

- `trait_list.txt`
- `traits/`
- `assoc/`
- `run_emmax_commands.sh`
- `plot_plan.txt`

## 需要时再读取

- 输入文件与列要求：`references/input_contract.md`
- 结果目录组织建议：`references/workflow.md`
