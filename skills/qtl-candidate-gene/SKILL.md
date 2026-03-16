---
name: qtl-candidate-gene
description: Use when GWAS summary statistics are available and you need to collapse significant SNPs into QTL regions, refine intervals with LD or fixed windows, annotate nearby genes, and export candidate gene tables for breeding interpretation; 适用于从 GWAS 结果出发完成 QTL 解释和候选基因提取。
---

# qtl-candidate-gene

## 何时使用

- 已经拿到 GWAS 显著位点或完整 summary statistics
- 需要从 SNP 层面上升到 QTL 区间层面
- 需要结合基因注释、LD 或固定窗口整理候选基因

## 快速开始

- 按 `references/input_contract.md` 准备 GWAS 结果和 gene annotation
- 复制并修改 `assets/qtl_candidate_gene.env.template`
- 运行 `bash scripts/run_qtl_candidate_gene.sh --config assets/qtl_candidate_gene.env.template`

## 默认工作流

1. 过滤出达到阈值的显著位点
2. 以 trait 和染色体为单位整理 lead SNP
3. 用 LD 表或固定窗口生成 QTL 区间
4. 将区间与 gene annotation 交叉，收集候选基因
5. 导出 QTL summary 与 candidate gene table

## 主要产物

- `significant_hits.tsv`
- `workflow_plan.txt`
- `qtl_regions/`
- `candidates/`

## 需要时再读取

- 输入结果表要求：`references/input_contract.md`
- QTL 与候选基因整理策略：`references/workflow.md`
