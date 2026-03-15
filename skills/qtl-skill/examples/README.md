# QTL Plotting Examples

This directory contains two plotting demonstrations for the QTL plotting layer.

## Example sets

### `basic/`

Contains the minimum tables needed to produce:

- a global Manhattan panel
- one local Manhattan panel
- one gene model panel without highlighted candidate genes

### `rich/`

Contains a richer example with:

- multiple loci
- LD-colored local points
- highlighted and labeled genes

## Files

- `chrom_sizes.tsv`
- `global_manhattan.tsv`
- `qtl_regions.tsv`
- `local_manhattan.tsv`
- `gene_models.tsv`
- `highlight_genes.tsv` (rich example only)
- rendered plot image(s)

## Renderers used

- `basic/basic_qtl_summary.png`
- `rich/rich_qtl_summary.png`

Both examples are generated through the default Python plotting workflow exposed by:

- `scripts/run_qtl_plot.sh`
