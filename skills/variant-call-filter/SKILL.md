---
name: variant-call-filter
description: Use when you already have aligned BAM files and need cohort-level variant calling, normalization, and filtering to obtain a clean VCF for downstream population genetics or GWAS; 适用于 BAM 已准备好、需要做变异检测与过滤清洗的场景。
---

# variant-call-filter

## 何时使用

- 上游已经完成 read QC 与 reference alignment
- 目标是得到高质量、可追溯的 cohort VCF
- 需要同时保留 raw call、normalized VCF 和 filtered VCF

## 快速开始

- 按 `references/input_contract.md` 准备 `bam_list.tsv`
- 复制并修改 `assets/variant_call_filter.env.template`
- 运行 `bash scripts/run_variant_call_filter.sh --config assets/variant_call_filter.env.template`

## 默认工作流

1. 校验 BAM 清单和参考基因组
2. 根据项目策略选择 `bcftools` 或 `GATK` caller
3. 对原始 VCF 做 left-normalization、拆分多等位位点和基本统计
4. 依据 `QUAL`、深度、缺失率和等位频率进行过滤
5. 输出过滤报告、位点统计和最终 cohort VCF

## 主要产物

- `bam_manifest.tsv`
- `bam_paths.list`
- `workflow_plan.txt`
- `raw_calls/`
- `normalized/`
- `filtered/`

## 需要时再读取

- BAM 清单格式：`references/input_contract.md`
- caller 与过滤阈值建议：`references/workflow.md`
