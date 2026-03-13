#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
data_dir="${GWAS_DATA_DIR:-${workspace_dir}/data}"
outdir="${RESQ_OUTDIR:-${workspace_dir}/analysis/resequencing_prep}"
logdir="${outdir}/logs"
mkdir -p "${outdir}" "${logdir}" "${outdir}/bam" "${outdir}/vcf" "${outdir}/genotype"

manifest="${RESQ_SAMPLE_MANIFEST:-${data_dir}/samples.tsv}"
reference="${RESQ_REFERENCE_FASTA:-${data_dir}/reference.fa}"
threads="${RESQ_THREADS:-4}"
num_pcs="${RESQ_NUM_PCS:-3}"
export RESQ_SAMPLE_MANIFEST="${manifest}"
export RESQ_OUTDIR="${outdir}"

aligner="${BWA_MEM2_BIN:-$(command -v bwa-mem2 2>/dev/null || true)}"
if [[ -z "${aligner}" ]]; then
  aligner="${BWA_BIN:-$(command -v bwa 2>/dev/null || true)}"
fi
samtools_bin="${SAMTOOLS_BIN:-$(command -v samtools 2>/dev/null || true)}"
bcftools_bin="${BCFTOOLS_BIN:-$(command -v bcftools 2>/dev/null || true)}"
plink_bin="${PLINK_BIN:-$(command -v plink 2>/dev/null || true)}"
emmax_kin_bin="${EMMAX_KIN_BIN:-$(command -v emmax-kin-intel64 2>/dev/null || true)}"

for tool_name in aligner samtools_bin bcftools_bin plink_bin emmax_kin_bin; do
  if [[ -z "${!tool_name}" ]]; then
    echo "Required tool for resequencing preparation is missing: ${tool_name}" >&2
    exit 1
  fi
done

if [[ ! -f "${reference}.fai" ]]; then
  "${samtools_bin}" faidx "${reference}"
fi
if [[ ! -f "${reference}.bwt" && ! -f "${reference}.0123" ]]; then
  "${aligner}" index "${reference}"
fi

python3 - <<'PY'
import pandas as pd
from pathlib import Path
import os
manifest = Path(os.environ["RESQ_SAMPLE_MANIFEST"])
df = pd.read_csv(manifest, sep="\t")
required = {"sample_id", "read1"}
if not required.issubset(df.columns):
    raise ValueError("samples.tsv must contain sample_id and read1 columns.")
PY

python3 - <<'PY'
import os
from pathlib import Path
import pandas as pd
manifest = pd.read_csv(os.environ["RESQ_SAMPLE_MANIFEST"], sep="\t")
Path(os.environ["RESQ_OUTDIR"]).mkdir(parents=True, exist_ok=True)
manifest["sample_id"].to_csv(Path(os.environ["RESQ_OUTDIR"]) / "sample_ids.txt", index=False, header=False)
PY

while IFS=$'\t' read -r sample_id read1 read2; do
  [[ "${sample_id}" == "sample_id" ]] && continue
  bam="${outdir}/bam/${sample_id}.sorted.bam"
  if [[ -f "${bam}" && -f "${bam}.bai" ]]; then
    continue
  fi
  if [[ -n "${read2:-}" ]]; then
    "${aligner}" mem -t "${threads}" "${reference}" "${read1}" "${read2}" | "${samtools_bin}" sort -@ "${threads}" -o "${bam}" -
  else
    "${aligner}" mem -t "${threads}" "${reference}" "${read1}" | "${samtools_bin}" sort -@ "${threads}" -o "${bam}" -
  fi
  "${samtools_bin}" index "${bam}"
done < "${manifest}"

find "${outdir}/bam" -name '*.sorted.bam' | sort > "${outdir}/bam/bam.list"

raw_vcf="${outdir}/vcf/raw_variants.vcf.gz"
snps_vcf="${outdir}/vcf/snps_only.vcf.gz"

if [[ ! -f "${raw_vcf}" ]]; then
  "${bcftools_bin}" mpileup -f "${reference}" -b "${outdir}/bam/bam.list" -Ou | \
    "${bcftools_bin}" call -mv -Oz -o "${raw_vcf}"
  "${bcftools_bin}" index "${raw_vcf}"
fi

if [[ ! -f "${snps_vcf}" ]]; then
  "${bcftools_bin}" view -m2 -M2 -v snps -Oz -o "${snps_vcf}" "${raw_vcf}"
  "${bcftools_bin}" index "${snps_vcf}"
fi

geno_prefix="${outdir}/genotype/genotype_panel"

if [[ ! -f "${geno_prefix}.bed" ]]; then
  "${plink_bin}" --vcf "${snps_vcf}" --double-id --allow-extra-chr --make-bed --out "${geno_prefix}"
fi

if [[ ! -f "${geno_prefix}.tped" ]]; then
  "${plink_bin}" --vcf "${snps_vcf}" --double-id --allow-extra-chr --recode transpose --out "${geno_prefix}"
fi

if [[ ! -f "${outdir}/genotype/kinship_matrix.tsv" ]]; then
  "${emmax_kin_bin}" "${geno_prefix}.tped" -o "${outdir}/genotype/kinship_matrix.tsv"
fi

if [[ ! -f "${outdir}/genotype/genotype_panel_prune.prune.in" ]]; then
  "${plink_bin}" --bfile "${geno_prefix}" --indep-pairwise 50 5 0.2 --allow-no-sex --out "${outdir}/genotype/genotype_panel_prune"
fi

if [[ ! -f "${outdir}/genotype/genotype_panel_pca.eigenvec" ]]; then
  "${plink_bin}" --bfile "${geno_prefix}" --extract "${outdir}/genotype/genotype_panel_prune.prune.in" --pca 10 --allow-no-sex --out "${outdir}/genotype/genotype_panel_pca"
fi

python3 "${script_dir}/make_covariates_from_pca.py" --eigenvec "${outdir}/genotype/genotype_panel_pca.eigenvec" --out "${outdir}/genotype/covariates.tsv" --num-pcs "${num_pcs}" --include-intercept
