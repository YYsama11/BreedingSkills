# Workflow Notes

## Interval definition strategy

- If LD has already been computed, prefer LD-supported boundaries for QTL intervals
- If LD is not available yet, use a fixed window as the initial interval
- Nearby peaks on the same trait and chromosome should be merged or resolved to one lead signal

## Candidate gene guidance

- Keep all genes inside each QTL interval
- Reserve separate columns for representative genes, nearest genes, and functional notes
- If annotation databases are available, append functional keywords in dedicated columns

## Delivery standard

- Every QTL can be traced back to its lead SNP and significance threshold
- The candidate gene table records the applied window or LD rule
- All outputs are ready for downstream plotting or report assembly
