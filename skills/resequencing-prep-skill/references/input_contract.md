# Input Contract

## Required files

- `samples.tsv`
- `reference.fa`

## Sample manifest format

Required columns:

- `sample_id`
- `read1`

Optional column:

- `read2`

Example:

```text
sample_id  read1                  read2
S1         /path/S1_R1.fastq.gz   /path/S1_R2.fastq.gz
S2         /path/S2_R1.fastq.gz   /path/S2_R2.fastq.gz
```

Single-end data is allowed if `read2` is missing or empty.

## Environment variables

- `GWAS_WORKSPACE_DIR`
- `GWAS_DATA_DIR`
- `RESQ_SAMPLE_MANIFEST`
- `RESQ_REFERENCE_FASTA`
- `RESQ_THREADS`
- `RESQ_NUM_PCS`
- `BWA_BIN` or `BWA_MEM2_BIN`
- `SAMTOOLS_BIN`
- `BCFTOOLS_BIN`
- `PLINK_BIN`
- `EMMAX_KIN_BIN`
