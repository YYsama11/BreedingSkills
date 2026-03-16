#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: bash scripts/run_gwas_visualization.sh --config path/to/gwas_run_visualization.env
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

required_vars=(TPED_PREFIX KINSHIP_FILE PHENO_TSV OUTDIR EMMAX_BIN)
for var_name in "${required_vars[@]}"; do
  [[ -n "${!var_name:-}" ]] || { echo "Missing variable: $var_name" >&2; exit 1; }
done

mkdir -p "$OUTDIR/traits" "$OUTDIR/assoc" "$OUTDIR/plots" "$OUTDIR/logs"

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

command_file="$OUTDIR/run_emmax_commands.sh"
printf '%s\n' '#!/usr/bin/env bash' 'set -euo pipefail' > "$command_file"

while IFS= read -r trait; do
  awk -F'\t' -v trait="$trait" '
  NR==1 {
    for (i = 1; i <= NF; i++) {
      idx[$i] = i
    }
    next
  }
  {
    print $1 "\t" $(idx[trait])
  }
  ' "$PHENO_TSV" > "$OUTDIR/traits/${trait}.pheno.tsv"

  if [[ -n "${COVARIATE_TSV:-}" ]]; then
    printf '%q ' "$EMMAX_BIN" -t "$TPED_PREFIX" -p "$OUTDIR/traits/${trait}.pheno.tsv" -k "$KINSHIP_FILE" -c "$COVARIATE_TSV" -o "$OUTDIR/assoc/${trait}" >> "$command_file"
  else
    printf '%q ' "$EMMAX_BIN" -t "$TPED_PREFIX" -p "$OUTDIR/traits/${trait}.pheno.tsv" -k "$KINSHIP_FILE" -o "$OUTDIR/assoc/${trait}" >> "$command_file"
  fi
  printf '\n' >> "$command_file"
done < "$OUTDIR/trait_list.txt"

chmod +x "$command_file"

cat > "$OUTDIR/plot_plan.txt" <<EOF
Expected per-trait summary columns:
trait	chrom	pos	snp_id	pvalue	beta	stderr

Expected visualization outputs:
plots/<trait>.manhattan.png
plots/<trait>.qq.png
assoc/<trait>.lead_snps.tsv

Reference plotting implementation:
skills/gwas-run-visualization/scripts/plot_manhattan_qq_reference.R

Reference style notes:
skills/gwas-run-visualization/references/plot_style.md
EOF

printf 'Prepared %s and %s\n' "$OUTDIR/trait_list.txt" "$command_file"
