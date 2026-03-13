# Observed Reference Run Summary

This reference summarizes the concrete large-scale run that motivated this skill.

## Input scale

- 112 target accessions
- 3 replicate measurements per accession
- 1461 molecular phenotypes
- 48 phenotype classes
- 4,033,947 SNPs in the genotype panel

## Generated phenotype layers

- molecule: 1461
- class: 48
- superclass: 6
- biological role: 5
- chain length bin: 4
- unsaturation bin: 4
- chain count bin: 5
- total abundance: 1

Total traits tested: **1534**

## EMMAX GWAS outcome

- 1534 EMMAX result sets produced
- 1534 Manhattan plots
- 1534 QQ plots
- 348 traits with at least one significant SNP
- 89,802 unique significant SNPs

## QTL / locus summary

- 552 global nonredundant lead loci after 500 kb collapsing
- 552 QTL regions
- 464 hotspot rows under the current hotspot rule
- 621 unique candidate genes
- 16 candidate genes matched lipid-related keyword rules after annotation backfill

## Population summary

- 3 inferred subgroups from genotype PCA
- SNP-level significant allele-frequency summary reduced to 89,802 rows after optimization

## Known limitations from this run

- Many loci collapse to lead-only intervals (`qtl_size_kb = 0`) under the chosen LD threshold.
- Non-synonymous SNP calling was not supported because the genotype alleles were encoded as `1/2`, not nucleotide bases.
- Hotspot counts are sensitive to very large shared-signal blocks and may need stricter user-defined thresholds for publication.
