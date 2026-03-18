<p align="center">
  <img src="BreedingSkills_logo_hd.png" alt="BreedingSkills logo" width="360">
</p>

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

## Agent 安装方式

BreedingSkills 不需要部署成在线服务；如果你想让 agent 使用它，核心是把这些 skill 安装到对应 agent 的本地 skills 目录中。

### Codex CLI

```bash
git clone https://github.com/YYsama11/BreedingSkills.git
cd BreedingSkills
bash install-codex.sh
```

### Claude Code

```bash
bash install-claude.sh
```

### Gemini CLI

```bash
bash install-gemini.sh
```

### OpenClaw

```bash
bash install-openclaw.sh
```

## 安装器能力

- `--list`：列出可安装 skill
- `--validate`：安装前校验 skill 结构
- `--skills read-qc-alignment,gwas-run-visualization`：只安装指定 skill
- `--project`：安装到当前项目目录，而不是全局目录
- `--update`：只更新源仓库里发生变化的 skill
- `--uninstall`：卸载此前安装过的 BreedingSkills
- `--dry-run`：仅预览安装结果

示例：

```bash
bash install-codex.sh --skills gwas-run-visualization,qtl-candidate-gene --project
```
