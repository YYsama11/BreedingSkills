# BreedingSkills

面向育种重测序、GWAS 与 QTL 解释流程的技能仓库。

当前版本先搭建 5 个核心 skill 的标准化骨架，所有 skill 都遵循同一目录约定：

- `SKILL.md`：核心说明文件
- `scripts/`：示例脚本
- `references/`：按需加载的参考资料
- `assets/`：模板文件

## Skills

- `skills/read-qc-alignment`：原始 FASTQ 质控、过滤与比对
- `skills/variant-call-filter`：变异检测、标准化与过滤
- `skills/genotype-kinship-prep`：EMMAX 输入、协变量与 kinship 准备
- `skills/gwas-run-visualization`：EMMAX 关联分析与结果整理
- `skills/qtl-candidate-gene`：QTL 解释与候选基因提取

## 入口

- 技能总览：`skills/INDEX.md`
- 每个 skill 的具体用法：对应目录下的 `SKILL.md`
