#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="${GWAS_WORKSPACE_DIR:-$(cd "${script_dir}/.." && pwd)}"
data_dir="${GWAS_DATA_DIR:-${workspace_dir}/data}"
outdir="${GWAS_GENOTYPE_OUTDIR:-${workspace_dir}/analysis/genotype}"
logdir="${GWAS_GENOTYPE_LOGDIR:-${workspace_dir}/analysis/logs/genotype}"
mkdir -p "${outdir}" "${logdir}"

plink_bin="${PLINK_BIN:-$(command -v plink 2>/dev/null || true)}"
if [[ -z "${plink_bin}" ]]; then
  echo "PLINK_BIN is not set and plink was not found in PATH." >&2
  exit 1
fi
threads="${GWAS_THREADS:-4}"
id_file="${GWAS_ID_FILE:-${data_dir}/sample_ids.txt}"
geno_prefix="${GWAS_GENO_PREFIX:-${data_dir}/genotype_panel}"
covariate_file="${GWAS_COVARIATE_FILE:-${data_dir}/covariates.tsv}"
kinship_file="${GWAS_KINSHIP_FILE:-${data_dir}/kinship_matrix.tsv}"

keep_file="${outdir}/keep_subset.txt"
awk '{print $1 "\t" $1}' "${id_file}" > "${keep_file}"

export GWAS_WORKSPACE_DIR="${workspace_dir}"
export GWAS_DATA_DIR="${data_dir}"
export GWAS_ID_FILE="${id_file}"
export GWAS_GENO_PREFIX="${geno_prefix}"
export GWAS_COVARIATE_FILE="${covariate_file}"
export GWAS_KINSHIP_FILE="${kinship_file}"
export GWAS_GENOTYPE_OUTDIR="${outdir}"

python3 - <<'PY'
import os
from pathlib import Path
import numpy as np
import pandas as pd

workspace_dir = Path(os.environ["GWAS_WORKSPACE_DIR"])
data_dir = Path(os.environ["GWAS_DATA_DIR"])
id_file = Path(os.environ["GWAS_ID_FILE"])
geno_prefix = os.environ["GWAS_GENO_PREFIX"]
covariate_file = Path(os.environ["GWAS_COVARIATE_FILE"])
kinship_file = Path(os.environ["GWAS_KINSHIP_FILE"])
outdir = Path(os.environ["GWAS_GENOTYPE_OUTDIR"])
outdir.mkdir(parents=True, exist_ok=True)

ids = [line.strip() for line in id_file.read_text().splitlines() if line.strip()]
tfam = pd.read_csv(f"{geno_prefix}.tfam", sep=r"\s+", header=None, names=["FID", "IID", "PID", "MID", "SEX", "PHENO"])
cov = pd.read_csv(covariate_file, sep=r"\s+", header=None)

subset_tfam = tfam[tfam["IID"].isin(ids)].copy()
subset_tfam = subset_tfam.set_index("IID").loc[ids].reset_index()
subset_tfam["FID"] = subset_tfam["IID"]
subset_tfam = subset_tfam[["FID", "IID", "PID", "MID", "SEX", "PHENO"]]
subset_tfam.to_csv(outdir / "sample_order.tsv", sep="\t", index=False)

cov.columns = ["FID", "IID"] + [f"COV{i}" for i in range(1, cov.shape[1] - 1)]
subset_cov = cov[cov["IID"].isin(ids)].copy().set_index("IID").loc[ids].reset_index()
subset_cov["FID"] = subset_cov["IID"]
subset_cov = subset_cov[["FID", "IID"] + [col for col in subset_cov.columns if col.startswith("COV")]]
subset_cov.to_csv(outdir / "covariates_subset.tsv", sep="\t", index=False)

kinship = np.loadtxt(kinship_file, dtype=float)
tfam_ids = tfam["IID"].tolist()
id_to_idx = {sample_id: idx for idx, sample_id in enumerate(tfam_ids)}
indices = [id_to_idx[sample_id] for sample_id in ids]
subset_kinship = kinship[np.ix_(indices, indices)]
np.savetxt(outdir / "kinship_subset.tsv", subset_kinship, fmt="%.10f", delimiter="\t")
PY

"${plink_bin}" \
  --tfile "${geno_prefix}" \
  --keep "${keep_file}" \
  --make-bed \
  --allow-no-sex \
  --threads "${threads}" \
  --out "${outdir}/target_subset" \
  > "${logdir}/plink_make_bed.log" 2>&1

"${plink_bin}" \
  --tfile "${geno_prefix}" \
  --keep "${keep_file}" \
  --recode transpose \
  --allow-no-sex \
  --threads "${threads}" \
  --out "${outdir}/target_subset_tped" \
  > "${logdir}/plink_make_tped.log" 2>&1

"${plink_bin}" \
  --bfile "${outdir}/target_subset" \
  --freq \
  --allow-no-sex \
  --threads "${threads}" \
  --out "${outdir}/target_subset_freq" \
  > "${logdir}/plink_freq.log" 2>&1

"${plink_bin}" \
  --bfile "${outdir}/target_subset" \
  --indep-pairwise 50 5 0.2 \
  --allow-no-sex \
  --threads "${threads}" \
  --out "${outdir}/target_subset_prune" \
  > "${logdir}/plink_prune.log" 2>&1

"${plink_bin}" \
  --bfile "${outdir}/target_subset" \
  --extract "${outdir}/target_subset_prune.prune.in" \
  --pca 10 \
  --allow-no-sex \
  --threads "${threads}" \
  --out "${outdir}/target_subset_pca" \
  > "${logdir}/plink_pca.log" 2>&1
