# BreedingSkills

<p align="center">
  <a href="#english">English</a> |
  <a href="#简体中文">简体中文</a>
</p>

<a id="english"></a>
<details open>
<summary><strong>English</strong></summary>

BreedingSkills is a skill repository for breeding-oriented resequencing, GWAS, and QTL interpretation workflows.

The current version provides five core skills, all following the same layout:

- `SKILL.md`: required core instructions
- `scripts/`: optional example scripts
- `references/`: optional reference material
- `assets/`: optional templates

## Skills

- `skills/read-qc-alignment`: raw FASTQ quality control, trimming, and read alignment
- `skills/variant-call-filter`: variant calling, normalization, and filtering
- `skills/genotype-kinship-prep`: EMMAX input preparation, covariates, and kinship
- `skills/gwas-run-visualization`: EMMAX GWAS execution and result organization
- `skills/qtl-candidate-gene`: QTL interpretation and candidate gene collection

## Entry Points

- Skill index: `skills/INDEX.md`
- Per-skill instructions: each skill folder contains its own `SKILL.md`

</details>

<a id="简体中文"></a>
<details>
<summary><strong>简体中文</strong></summary>

BreedingSkills 是一个面向育种重测序、GWAS 和 QTL 解析流程的技能仓库。

当前版本包含 5 个核心 skill，且都遵循统一目录结构：

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

</details>
