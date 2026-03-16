---
name: genotype-kinship-prep
description: Use when you need to transform filtered genotype data plus phenotype and covariate tables into EMMAX-ready inputs including TPED/TFAM, kinship, PCA covariates, and sample-intersected trait matrices; 适用于为 emmax-intel64 准备 GWAS 输入文件的场景。
---

# genotype-kinship-prep

## 何时使用

- 已经拥有过滤后的 VCF 或 PLINK 数据
- 需要整理样本交集、协变量、PCA 与 kinship
- 下游关联分析软件固定为 `emmax-intel64`

## 快速开始

- 按 `references/input_contract.md` 准备基因型和表型矩阵
- 复制并修改 `assets/genotype_kinship_prep.env.template`
- 运行 `bash scripts/run_genotype_kinship_prep.sh --config assets/genotype_kinship_prep.env.template`

## 默认工作流

1. 校验 genotype、phenotype 和 sample ID 是否一致
2. 将 VCF 转换为 EMMAX 需要的 `tped/tfam`
3. 构建 PLINK 二进制文件并计算 PCA 协变量
4. 生成 kinship 矩阵
5. 提取 trait 列表，形成批量运行的 trait manifest

## 主要产物

- `trait_list.txt`
- `workflow_plan.txt`
- `emmax/`
- `covariates/`
- `intermediate/`

## 需要时再读取

- 输入矩阵格式：`references/input_contract.md`
- `emmax-intel64` 准备细节：`references/workflow.md`
