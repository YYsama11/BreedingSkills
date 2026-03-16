# Skill Index

## 统一结构

- `SKILL.md`：定义触发条件、工作流、输入输出和引用资源
- `scripts/`：尽量放确定性高、可复用的启动脚本
- `references/`：把输入约束、流程细节放到按需读取的文档中
- `assets/`：配置模板、表头模板或其他可复制资源

## 当前技能

- `read-qc-alignment`：从原始 FASTQ 到排序 BAM 与比对质控
- `variant-call-filter`：从 BAM 到标准化、过滤后的 cohort VCF
- `genotype-kinship-prep`：把 VCF/表型整理为 `emmax-intel64` 可直接使用的输入
- `gwas-run-visualization`：批量运行 EMMAX 并整理 Manhattan / QQ / lead SNP 结果
- `qtl-candidate-gene`：把显著位点整理为 QTL 区间并收集候选基因
