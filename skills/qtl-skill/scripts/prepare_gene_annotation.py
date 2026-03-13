#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


ATTR_RE = re.compile(r"([^=;]+)=([^;]+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create BED-like annotation resources from GFF.")
    parser.add_argument("--gff", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--promoter-length", type=int, default=2000)
    return parser.parse_args()


def parse_attr(raw: str) -> dict[str, str]:
    return {key: value for key, value in ATTR_RE.findall(raw)}


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = []
    mrna_to_gene = {}
    gene_to_mrna_annot: dict[str, dict[str, str]] = {}
    exons = []
    cds = []

    with open(args.gff) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, source, feature, start, end, score, strand, phase, attrs_raw = line.rstrip("\n").split("\t")
            attrs = parse_attr(attrs_raw)
            start_i = int(start)
            end_i = int(end)
            if feature == "gene":
                gene_id = attrs.get("ID")
                genes.append(
                    {
                        "chrom": chrom,
                        "start": start_i,
                        "end": end_i,
                        "gene_id": gene_id,
                        "strand": strand,
                        "source": source,
                        "nr_annotation": attrs.get("nr_annotation", ""),
                        "pfam": attrs.get("pfam", ""),
                    }
                )
            elif feature == "mRNA":
                mrna_id = attrs.get("ID")
                parent = attrs.get("Parent")
                if mrna_id and parent:
                    mrna_to_gene[mrna_id] = parent
                    if parent not in gene_to_mrna_annot:
                        gene_to_mrna_annot[parent] = {
                            "nr_annotation": attrs.get("nr_annotation", ""),
                            "pfam": attrs.get("pfam", ""),
                            "source": source,
                        }
            elif feature in {"exon", "five_prime_UTR", "three_prime_UTR"}:
                parent = attrs.get("Parent")
                exons.append((chrom, start_i, end_i, strand, parent))
            elif feature == "CDS":
                parent = attrs.get("Parent")
                cds.append((chrom, start_i, end_i, strand, parent))

    gene_df = pd.DataFrame(genes).dropna(subset=["gene_id"]).drop_duplicates(subset=["gene_id"]).copy()
    if not gene_df.empty:
        gene_df["nr_annotation"] = gene_df.apply(
            lambda row: row["nr_annotation"] if row["nr_annotation"] else gene_to_mrna_annot.get(row["gene_id"], {}).get("nr_annotation", ""),
            axis=1,
        )
        gene_df["pfam"] = gene_df.apply(
            lambda row: row["pfam"] if row["pfam"] else gene_to_mrna_annot.get(row["gene_id"], {}).get("pfam", ""),
            axis=1,
        )
    gene_df["bed_start"] = gene_df["start"] - 1
    gene_df["promoter_start"] = gene_df.apply(
        lambda row: max(0, row["start"] - args.promoter_length - 1) if row["strand"] == "+" else row["end"],
        axis=1,
    )
    gene_df["promoter_end"] = gene_df.apply(
        lambda row: row["start"] - 1 if row["strand"] == "+" else row["end"] + args.promoter_length,
        axis=1,
    )
    gene_df["downstream_start"] = gene_df.apply(
        lambda row: row["end"] if row["strand"] == "+" else max(0, row["start"] - args.promoter_length - 1),
        axis=1,
    )
    gene_df["downstream_end"] = gene_df.apply(
        lambda row: row["end"] + args.promoter_length if row["strand"] == "+" else row["start"] - 1,
        axis=1,
    )

    gene_df[["chrom", "bed_start", "end", "gene_id", "strand"]].to_csv(
        outdir / "genes.bed", sep="\t", header=False, index=False
    )
    gene_df[["chrom", "promoter_start", "promoter_end", "gene_id", "strand"]].to_csv(
        outdir / "promoters_2kb.bed", sep="\t", header=False, index=False
    )
    gene_df[["chrom", "downstream_start", "downstream_end", "gene_id", "strand"]].to_csv(
        outdir / "downstream_2kb.bed", sep="\t", header=False, index=False
    )
    gene_df.to_csv(outdir / "gene_metadata.tsv", sep="\t", index=False)

    exon_records = []
    for chrom, start_i, end_i, strand, parent in exons:
        gene_id = mrna_to_gene.get(parent, parent)
        exon_records.append((chrom, start_i - 1, end_i, gene_id, strand))
    if exon_records:
        pd.DataFrame(exon_records).drop_duplicates().to_csv(
            outdir / "exons.bed", sep="\t", header=False, index=False
        )

    cds_records = []
    for chrom, start_i, end_i, strand, parent in cds:
        gene_id = mrna_to_gene.get(parent, parent)
        cds_records.append((chrom, start_i - 1, end_i, gene_id, strand))
    if cds_records:
        pd.DataFrame(cds_records).drop_duplicates().to_csv(
            outdir / "cds.bed", sep="\t", header=False, index=False
        )


if __name__ == "__main__":
    main()
