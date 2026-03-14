# Plot Input Contract

## Minimal plotting tables

### `chrom_sizes`

- `chrom`
- `length_bp`

Controls the chromosome layout in the top genome-wide Manhattan panel.

### `global_manhattan`

- `snp_id`
- `chrom`
- `bp`
- `p_value` or `mlog10_p`

Provides the scatter points for the top genome-wide Manhattan panel.

### `qtl_regions`

- `locus_id`
- `chrom`
- `lead_snp_id`
- `lead_bp`
- `lead_p`
- `qtl_start`
- `qtl_end`
- `panel_start`
- `panel_end`
- `qtl_type`

Defines the local zoom windows and lead SNP metadata for each locus.

### `local_manhattan`

- `locus_id`
- `snp_id`
- `chrom`
- `bp`
- `mlog10_p` or `p_value`
- `r2` (optional)
- `is_lead` (optional)

Provides the scatter points for each local Manhattan panel.

### `gene_models`

- `gene_id`
- `chrom`
- `start_bp`
- `end_bp`
- `strand`

Provides gene structures for the lower panel of each locus.

### `highlight_genes` (optional)

- `locus_id`
- `gene_id`
- `rank` (optional)
- `label` (optional)

If missing, all genes are drawn in gray without highlighted labels.

## Rendering behavior

- If `chrom_sizes` is missing, chromosome sizes are inferred from `global_manhattan` and `qtl_regions`.
- If `r2` is missing, local points are drawn in neutral gray.
- If `gene_models` is missing, the lower panel is drawn as an empty annotation track with a message.
