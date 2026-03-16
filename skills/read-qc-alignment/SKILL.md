---
name: read-qc-alignment
description: Use when starting from raw FASTQ reads and you need resequencing QC, trimming, reference alignment, sorted BAM generation, and mapping QC for breeding projects; 适用于从原始 FASTQ 开始完成质控、过滤、比对和 BAM 质控。
---

# read-qc-alignment

## 何时使用

- 输入仍然是原始测序数据，通常为双端 `FASTQ.gz`
- 目标是生成后续变异检测可用的高质量 BAM
- 需要在流程里同时保留原始质控和过滤后质控结果

## 快速开始

- 按 `references/input_contract.md` 准备样本表与参考基因组
- 复制并修改 `assets/read_qc_alignment.env.template`
- 运行 `bash scripts/run_read_qc_alignment.sh --config assets/read_qc_alignment.env.template`

## 默认工作流

1. 校验样本表、FASTQ 路径和参考基因组索引
2. 运行原始读段质控并汇总 FastQC 指标
3. 使用 `fastp` 或等价工具进行接头和低质量碱基过滤
4. 使用 `bwa mem` 或兼容短读长比对器完成参考基因组比对
5. 对 BAM 进行排序、索引和基础 mapping QC
6. 输出可供下游 variant calling 使用的 BAM 清单

## 主要产物

- `sample_manifest.tsv`
- `workflow_plan.txt`
- `qc_raw/` 与 `qc_trimmed/`
- `bam/` 与 `stats/`

## 需要时再读取

- 输入列和目录约束：`references/input_contract.md`
- 工具选择与关键检查点：`references/workflow.md`
