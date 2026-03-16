#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


GENOMEWIDE_P_DEFAULT = 0.05 / 4033947


@dataclass(frozen=True)
class Gene:
    gene_id: str
    chrom: int
    start: int
    end: int
    strand: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build fixed-window QTL summaries and plotting tables.")
    parser.add_argument("--gwas-table", required=True)
    parser.add_argument("--gene-annotation", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--p-threshold", type=float, required=True)
    parser.add_argument("--window-bp", type=int, required=True)
    parser.add_argument("--trait-label", default="")
    parser.add_argument("--reference-fai", default="")
    parser.add_argument("--max-loci", type=int, default=3)
    parser.add_argument("--genomewide-p", type=float, default=GENOMEWIDE_P_DEFAULT)
    return parser.parse_args()


def chrom_num(raw: str) -> int:
    match = re.search(r"(\d+)$", str(raw).strip())
    if not match:
        raise ValueError(f"Cannot parse chromosome number from: {raw}")
    return int(match.group(1))


def read_gff_genes(path: Path) -> tuple[dict[int, list[Gene]], dict[int, list[int]]]:
    genes_by_chr: dict[int, list[Gene]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            try:
                chrom = chrom_num(fields[0])
            except ValueError:
                continue
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attrs = {}
            for part in fields[8].split(";"):
                if "=" in part:
                    key, value = part.split("=", 1)
                    attrs[key] = value
            gene_id = attrs.get("ID", f"gene_{chrom}_{start}_{end}")
            genes_by_chr[chrom].append(Gene(gene_id=gene_id, chrom=chrom, start=start, end=end, strand=strand))
    starts_by_chr: dict[int, list[int]] = {}
    for chrom, genes in genes_by_chr.items():
        genes.sort(key=lambda item: (item.start, item.end, item.gene_id))
        starts_by_chr[chrom] = [gene.start for gene in genes]
    return genes_by_chr, starts_by_chr


def genes_in_region(genes_by_chr: dict[int, list[Gene]], starts_by_chr: dict[int, list[int]], chrom: int, start: int, end: int) -> list[Gene]:
    genes = genes_by_chr.get(chrom, [])
    starts = starts_by_chr.get(chrom, [])
    if not genes:
        return []
    index = bisect_left(starts, start)
    while index > 0 and genes[index - 1].end >= start:
        index -= 1
    out = []
    while index < len(genes) and genes[index].start <= end:
        gene = genes[index]
        if gene.end >= start:
            out.append(gene)
        index += 1
    return out


def nearest_genes(genes_by_chr: dict[int, list[Gene]], starts_by_chr: dict[int, list[int]], chrom: int, pos: int, limit: int = 3) -> list[Gene]:
    genes = genes_by_chr.get(chrom, [])
    starts = starts_by_chr.get(chrom, [])
    if not genes:
        return []
    index = bisect_left(starts, pos)
    candidates: list[tuple[int, Gene]] = []
    for j in range(max(0, index - 10), min(len(genes), index + 10)):
        gene = genes[j]
        distance = 0 if gene.start <= pos <= gene.end else min(abs(pos - gene.start), abs(pos - gene.end))
        candidates.append((distance, gene))
    candidates.sort(key=lambda item: (item[0], item[1].start, item[1].gene_id))
    return [gene for _, gene in candidates[:limit]]


def load_chrom_sizes(path: str) -> dict[int, int]:
    if not path:
        return {}
    sizes: dict[int, int] = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            sizes[chrom_num(fields[0])] = int(fields[1])
    return sizes


def detect_schema(fieldnames: list[str]) -> str:
    names = set(fieldnames)
    if {"trait", "chrom", "pos", "snp_id", "pvalue"}.issubset(names):
        return "gwas_summary"
    if {"SNP", "CHR", "BP", "P"}.issubset(names):
        return "pmap"
    raise ValueError(f"Unsupported GWAS table schema: {fieldnames}")


def standardize_row(row: dict[str, str], schema: str, trait_label: str) -> dict[str, object]:
    if schema == "gwas_summary":
        return {
            "trait": row["trait"] or trait_label,
            "chrom": chrom_num(row["chrom"]),
            "bp": int(float(row["pos"])),
            "snp_id": row["snp_id"],
            "pvalue": float(row["pvalue"]),
            "beta": row.get("beta", ""),
            "stderr": row.get("stderr", ""),
        }
    return {
        "trait": trait_label,
        "chrom": chrom_num(row["CHR"]),
        "bp": int(float(row["BP"])),
        "snp_id": row["SNP"],
        "pvalue": float(row["P"]),
        "beta": "",
        "stderr": "",
    }


def build_loci(sig_hits: list[dict[str, object]], cluster_bp: int, genomewide_p: float, max_loci: int) -> list[dict[str, object]]:
    loci: list[dict[str, object]] = []
    by_chrom: dict[int, list[dict[str, object]]] = defaultdict(list)
    for hit in sig_hits:
        by_chrom[int(hit["chrom"])].append(hit)
    for chrom, hits in by_chrom.items():
        work = sorted(hits, key=lambda item: (float(item["pvalue"]), int(item["bp"]), str(item["snp_id"])))
        while work:
            lead = work[0]
            lo = max(1, int(lead["bp"]) - cluster_bp)
            hi = int(lead["bp"]) + cluster_bp
            members = [row for row in work if lo <= int(row["bp"]) <= hi]
            work = [row for row in work if not (lo <= int(row["bp"]) <= hi)]
            loci.append(
                {
                    "chrom": chrom,
                    "lead_snp": str(lead["snp_id"]),
                    "lead_bp": int(lead["bp"]),
                    "lead_p": float(lead["pvalue"]),
                    "coarse_start": min(int(row["bp"]) for row in members),
                    "coarse_end": max(int(row["bp"]) for row in members),
                    "member_count": len(members),
                    "gw_member_count": sum(1 for row in members if float(row["pvalue"]) <= genomewide_p),
                    "members": members,
                }
            )
    loci.sort(key=lambda item: (float(item["lead_p"]), int(item["chrom"]), int(item["lead_bp"])))
    for index, locus in enumerate(loci[:max_loci], start=1):
        locus["locus_id"] = f"L{index:03d}"
        locus["qtl_type"] = "fixed_window"
        locus["qtl_start"] = max(1, int(locus["coarse_start"]) - cluster_bp)
        locus["qtl_end"] = int(locus["coarse_end"]) + cluster_bp
        locus["panel_start"] = max(1, int(locus["qtl_start"]) - cluster_bp)
        locus["panel_end"] = int(locus["qtl_end"]) + cluster_bp
    return loci[:max_loci]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    gwas_path = Path(args.gwas_table)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    trait_label = args.trait_label or gwas_path.name.replace(".gwas_summary.tsv", "").replace(".all.QQ.pmap.txt", "")
    genes_by_chr, starts_by_chr = read_gff_genes(Path(args.gene_annotation))
    chrom_sizes = load_chrom_sizes(args.reference_fai)

    significant_hits: list[dict[str, object]] = []
    global_rows: list[dict[str, object]] = []

    with gwas_path.open("r", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit(f"Empty GWAS table: {gwas_path}")
        schema = detect_schema(reader.fieldnames)
        for row in reader:
            std = standardize_row(row, schema=schema, trait_label=trait_label)
            chrom = int(std["chrom"])
            bp = int(std["bp"])
            pvalue = max(float(std["pvalue"]), float.fromhex("0x1.0p-1022"))
            chrom_sizes[chrom] = max(chrom_sizes.get(chrom, 0), bp)
            global_rows.append(
                {
                    "snp_id": std["snp_id"],
                    "chrom": chrom,
                    "bp": bp,
                    "p_value": pvalue,
                    "mlog10_p": -math.log10(pvalue),
                }
            )
            if pvalue <= args.p_threshold:
                significant_hits.append(std)

    if not global_rows:
        raise SystemExit(f"No rows found in GWAS table: {gwas_path}")

    loci = build_loci(significant_hits, args.window_bp, args.genomewide_p, args.max_loci)
    if not loci:
        raise SystemExit(f"No significant loci found at p <= {args.p_threshold}")

    plot_dir = outdir / "plot_tables"
    qtl_rows: list[dict[str, object]] = []
    candidate_rows: list[dict[str, object]] = []
    qtl_plot_rows: list[dict[str, object]] = []
    highlight_rows: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    added_genes: set[tuple[str, str]] = set()

    for locus in loci:
        genes = genes_in_region(genes_by_chr, starts_by_chr, int(locus["chrom"]), int(locus["qtl_start"]), int(locus["qtl_end"]))
        if not genes:
            genes = nearest_genes(genes_by_chr, starts_by_chr, int(locus["chrom"]), int(locus["lead_bp"]), limit=10)
        genes = genes[:10]
        top_gene_ids = [gene.gene_id for gene in genes[:3]]

        qtl_rows.append(
            {
                "trait": trait_label,
                "locus_id": locus["locus_id"],
                "chrom": locus["chrom"],
                "lead_snp": locus["lead_snp"],
                "lead_bp": locus["lead_bp"],
                "lead_p": locus["lead_p"],
                "qtl_type": locus["qtl_type"],
                "coarse_start": locus["coarse_start"],
                "coarse_end": locus["coarse_end"],
                "qtl_start": locus["qtl_start"],
                "qtl_end": locus["qtl_end"],
                "qtl_span_bp": int(locus["qtl_end"]) - int(locus["qtl_start"]) + 1,
                "member_count": locus["member_count"],
                "gw_member_count": locus["gw_member_count"],
                "top_candidate_genes": "; ".join(top_gene_ids),
            }
        )
        qtl_plot_rows.append(
            {
                "locus_id": locus["locus_id"],
                "chrom": locus["chrom"],
                "lead_snp_id": locus["lead_snp"],
                "lead_bp": locus["lead_bp"],
                "lead_p": locus["lead_p"],
                "qtl_start": locus["qtl_start"],
                "qtl_end": locus["qtl_end"],
                "panel_start": locus["panel_start"],
                "panel_end": locus["panel_end"],
                "qtl_type": locus["qtl_type"],
            }
        )

        for rank, gene in enumerate(genes, start=1):
            distance = 0 if gene.start <= int(locus["lead_bp"]) <= gene.end else min(abs(int(locus["lead_bp"]) - gene.start), abs(int(locus["lead_bp"]) - gene.end))
            candidate_rows.append(
                {
                    "trait": trait_label,
                    "locus_id": locus["locus_id"],
                    "rank": rank,
                    "gene_id": gene.gene_id,
                    "chrom": gene.chrom,
                    "gene_start": gene.start,
                    "gene_end": gene.end,
                    "lead_snp": locus["lead_snp"],
                    "lead_bp": locus["lead_bp"],
                    "lead_p": locus["lead_p"],
                    "distance_bp": distance,
                    "in_qtl": int(gene.start <= int(locus["qtl_end"]) and gene.end >= int(locus["qtl_start"])),
                }
            )
            if rank <= 3:
                highlight_rows.append(
                    {
                        "locus_id": locus["locus_id"],
                        "gene_id": gene.gene_id,
                        "rank": rank,
                        "label": gene.gene_id,
                    }
                )
            key = (locus["locus_id"], gene.gene_id)
            if key not in added_genes:
                gene_rows.append(
                    {
                        "gene_id": gene.gene_id,
                        "chrom": gene.chrom,
                        "start_bp": gene.start,
                        "end_bp": gene.end,
                        "strand": gene.strand,
                    }
                )
                added_genes.add(key)

    local_rows: list[dict[str, object]] = []
    loci_by_chrom: dict[int, list[dict[str, object]]] = defaultdict(list)
    for locus in loci:
        loci_by_chrom[int(locus["chrom"])].append(locus)

    for row in global_rows:
        chrom = int(row["chrom"])
        bp = int(row["bp"])
        for locus in loci_by_chrom.get(chrom, []):
            if int(locus["panel_start"]) <= bp <= int(locus["panel_end"]):
                local_rows.append(
                    {
                        "locus_id": locus["locus_id"],
                        "snp_id": row["snp_id"],
                        "chrom": chrom,
                        "bp": bp,
                        "p_value": row["p_value"],
                        "mlog10_p": row["mlog10_p"],
                    }
                )

    write_tsv(
        outdir / f"{trait_label}.qtl_summary.tsv",
        ["trait", "locus_id", "chrom", "lead_snp", "lead_bp", "lead_p", "qtl_type", "coarse_start", "coarse_end", "qtl_start", "qtl_end", "qtl_span_bp", "member_count", "gw_member_count", "top_candidate_genes"],
        qtl_rows,
    )
    write_tsv(
        outdir / f"{trait_label}.top10_candidate_genes.tsv",
        ["trait", "locus_id", "rank", "gene_id", "chrom", "gene_start", "gene_end", "lead_snp", "lead_bp", "lead_p", "distance_bp", "in_qtl"],
        candidate_rows,
    )
    write_tsv(plot_dir / "chrom_sizes.tsv", ["chrom", "length_bp"], [{"chrom": chrom, "length_bp": chrom_sizes[chrom]} for chrom in sorted(chrom_sizes)])
    write_tsv(plot_dir / "global_manhattan.tsv", ["snp_id", "chrom", "bp", "p_value", "mlog10_p"], global_rows)
    write_tsv(plot_dir / "qtl_regions.tsv", ["locus_id", "chrom", "lead_snp_id", "lead_bp", "lead_p", "qtl_start", "qtl_end", "panel_start", "panel_end", "qtl_type"], qtl_plot_rows)
    write_tsv(plot_dir / "local_manhattan.tsv", ["locus_id", "snp_id", "chrom", "bp", "p_value", "mlog10_p"], local_rows)
    write_tsv(plot_dir / "gene_models.tsv", ["gene_id", "chrom", "start_bp", "end_bp", "strand"], gene_rows)
    write_tsv(plot_dir / "highlight_genes.tsv", ["locus_id", "gene_id", "rank", "label"], highlight_rows)
    (outdir / f"{trait_label}.status.json").write_text(json.dumps({"trait": trait_label, "selected_loci": len(loci), "candidate_genes": len(candidate_rows)}, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
