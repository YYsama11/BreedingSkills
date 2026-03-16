# BreedingSkills

[English](README.md) | [简体中文](README.zh-CN.md)

BreedingSkills 是一个面向育种重测序、GWAS 和 QTL 解析流程的技能仓库。

当前版本包含 5 个核心 skill，且每个 skill 都遵循统一目录结构：

- `SKILL.md`：必需的核心说明文件
- `scripts/`：可选的示例脚本
- `references/`：可选的参考资料
- `assets/`：可选的模板文件

## Skills

- `skills/read-qc-alignment`：原始 FASTQ 质控、过滤与比对
- `skills/variant-call-filter`：变异检测、标准化与过滤
- `skills/genotype-kinship-prep`：EMMAX 输入、协变量与 kinship 准备
- `skills/gwas-run-visualization`：EMMAX 关联分析与结果整理
- `skills/qtl-candidate-gene`：QTL 解释与候选基因提取

## 入口

- 技能索引：`skills/INDEX.md`
- 每个技能的说明：对应目录下的 `SKILL.md`
