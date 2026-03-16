#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/run_genotype_kinship_prep.sh --config path/to/genotype_kinship_prep.env
EOF
}

config=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      config="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

[[ -n "$config" ]] || { usage; exit 1; }
[[ -f "$config" ]] || { echo "Config not found: $config" >&2; exit 1; }

set -a
source "$config"
set +a

required_vars=(VCF_FILE PHENO_TSV OUTDIR PREFIX THREADS EMMAX_BIN KINSHIP_BIN PCA_COMPONENTS)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

mkdir -p "$OUTDIR/intermediate" "$OUTDIR/covariates" "$OUTDIR/emmax" "$OUTDIR/logs"

awk -F'\t' '
NR==1 {
  if ($1 != "sample_id") {
    print "phenotype header must start with sample_id" > "/dev/stderr"
    exit 2
  }
  for (i = 2; i <= NF; i++) {
    print $i
  }
}
' "$PHENO_TSV" > "$OUTDIR/trait_list.txt"

cat > "$OUTDIR/workflow_plan.txt" <<EOF
plink --vcf ${VCF_FILE} --make-bed --out ${OUTDIR}/intermediate/${PREFIX}
plink --vcf ${VCF_FILE} --recode12 --transpose --output-missing-genotype 0 --out ${OUTDIR}/emmax/${PREFIX}
plink --bfile ${OUTDIR}/intermediate/${PREFIX} --pca ${PCA_COMPONENTS} header tabs --out ${OUTDIR}/covariates/${PREFIX}
${KINSHIP_BIN} -v -s -d 10 ${OUTDIR}/emmax/${PREFIX}
${EMMAX_BIN} -h
EOF

printf 'Prepared %s and %s\n' "$OUTDIR/trait_list.txt" "$OUTDIR/workflow_plan.txt"
