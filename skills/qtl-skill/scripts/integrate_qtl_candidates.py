#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Integrate LD-derived QTL regions and candidate gene analyses.")
    parser.add_argument("--lead-snps", required=True)
    parser.add_argument("--ld-dir", required=True)
    parser.add_argument("--master-significant", required=True)
    parser.add_argument("--trait-summaries", required=False)
    parser.add_argument("--annotation-dir", required=False)
    parser.add_argument("--gene-annotation-file", required=False)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--module-membership", required=False)
    return parser.parse_args()


def normalize_chrom_name(name: str) -> str:
    text = str(name).strip()
    lowered = text.lower()
    if lowered.startswith("chr"):
        lowered = lowered[3:]
    if lowered.isdigit():
        lowered = str(int(lowered))
    return lowered


def load_ld_region(ld_file: Path, lead_row: pd.Series) -> dict[str, object]:
    if not ld_file.exists() or ld_file.stat().st_size == 0:
        return {
            "lead_snp": lead_row["snp_id"],
            "chrom": str(lead_row["chrom"]),
            "lead_pos": int(lead_row["pos"]),
            "qtl_start": int(lead_row["pos"]),
            "qtl_end": int(lead_row["pos"]),
            "qtl_size_bp": 0,
            "ld_partner_count": 1,
            "max_r2": 1.0,
        }
    ld = pd.read_csv(ld_file, sep=r"\s+")
    if ld.empty:
        return {
            "lead_snp": lead_row["snp_id"],
            "chrom": str(lead_row["chrom"]),
            "lead_pos": int(lead_row["pos"]),
            "qtl_start": int(lead_row["pos"]),
            "qtl_end": int(lead_row["pos"]),
            "qtl_size_bp": 0,
            "ld_partner_count": 1,
            "max_r2": 1.0,
        }
    high_ld = ld[ld["R2"] >= 0.6].copy()
    if high_ld.empty:
        high_ld = ld.copy()
    qtl_start = int(high_ld["BP_B"].min())
    qtl_end = int(high_ld["BP_B"].max())
    return {
        "lead_snp": lead_row["snp_id"],
        "chrom": str(lead_row["chrom"]),
        "lead_pos": int(lead_row["pos"]),
        "qtl_start": qtl_start,
        "qtl_end": qtl_end,
        "qtl_size_bp": qtl_end - qtl_start,
        "ld_partner_count": int(high_ld.shape[0]),
        "max_r2": float(high_ld["R2"].max()),
    }


def annotate_snps(sig: pd.DataFrame, genes: pd.DataFrame, promoters: pd.DataFrame, downstream: pd.DataFrame, exons: pd.DataFrame, cds: pd.DataFrame) -> pd.DataFrame:
    def build_lookup(df: pd.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        lookup = {}
        for chrom, sub in df.groupby("chrom"):
            ordered = sub.sort_values("start")
            lookup[chrom] = (
                ordered["start"].to_numpy(dtype=np.int64),
                ordered["end"].to_numpy(dtype=np.int64),
                ordered["gene_id"].astype(str).to_numpy(),
            )
        return lookup

    def query_interval(positions: np.ndarray, lookup_item: tuple[np.ndarray, np.ndarray, np.ndarray] | None) -> tuple[np.ndarray, np.ndarray]:
        matched = np.zeros(len(positions), dtype=bool)
        genes_out = np.array([""] * len(positions), dtype=object)
        if lookup_item is None:
            return matched, genes_out
        starts, ends, gene_ids = lookup_item
        idx = np.searchsorted(starts, positions, side="right") - 1
        valid = idx >= 0
        if not valid.any():
            return matched, genes_out
        idx_valid = idx[valid]
        pos_valid = positions[valid]
        within = pos_valid <= ends[idx_valid]
        if within.any():
            matched_idx = np.where(valid)[0][within]
            matched[matched_idx] = True
            genes_out[matched_idx] = gene_ids[idx_valid[within]]
        return matched, genes_out

    cds_lookup = build_lookup(cds)
    exon_lookup = build_lookup(exons)
    prom_lookup = build_lookup(promoters)
    downstream_lookup = build_lookup(downstream)
    gene_lookup = build_lookup(genes)

    parts = []
    for chrom_value, sub in sig.groupby("chrom"):
        chrom = str(chrom_value)
        positions = sub["pos"].to_numpy(dtype=np.int64)
        region = np.array(["intergenic"] * len(sub), dtype=object)
        gene_ids = np.array([""] * len(sub), dtype=object)

        matched, genes_hit = query_interval(positions, cds_lookup.get(chrom))
        region[matched] = "CDS"
        gene_ids[matched] = genes_hit[matched]

        remaining = region == "intergenic"
        if remaining.any():
            matched, genes_hit = query_interval(positions[remaining], exon_lookup.get(chrom))
            region[np.where(remaining)[0][matched]] = "exon"
            gene_ids[np.where(remaining)[0][matched]] = genes_hit[matched]

        remaining = region == "intergenic"
        if remaining.any():
            matched, genes_hit = query_interval(positions[remaining], prom_lookup.get(chrom))
            region[np.where(remaining)[0][matched]] = "promoter"
            gene_ids[np.where(remaining)[0][matched]] = genes_hit[matched]

        remaining = region == "intergenic"
        if remaining.any():
            matched, genes_hit = query_interval(positions[remaining], downstream_lookup.get(chrom))
            region[np.where(remaining)[0][matched]] = "downstream_2kb"
            gene_ids[np.where(remaining)[0][matched]] = genes_hit[matched]

        remaining = region == "intergenic"
        if remaining.any():
            matched, genes_hit = query_interval(positions[remaining], gene_lookup.get(chrom))
            region[np.where(remaining)[0][matched]] = "intron"
            gene_ids[np.where(remaining)[0][matched]] = genes_hit[matched]

        sub_out = sub.loc[:, ["snp_id", "trait_id", "chrom", "pos"]].copy()
        sub_out["region_type"] = region
        sub_out["gene_id"] = gene_ids
        parts.append(sub_out)

    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["snp_id", "trait_id", "chrom", "pos", "region_type", "gene_id"])


def build_hotspots(qtl: pd.DataFrame) -> pd.DataFrame:
    if qtl.empty:
        return pd.DataFrame(columns=["hotspot_id", "chrom", "start", "end", "trait_count", "trait_ids"])

    def extract_traits(row: pd.Series) -> list[str]:
        if "trait_id" in row and pd.notna(row["trait_id"]):
            return [str(row["trait_id"])]
        if "trait_ids" in row and pd.notna(row["trait_ids"]):
            return [item for item in str(row["trait_ids"]).split(",") if item]
        return []

    hotspot_rows = []
    hotspot_id = 1
    for chrom, sub in qtl.sort_values(["chrom", "qtl_start", "qtl_end"]).groupby("chrom"):
        current_start = None
        current_end = None
        current_traits: list[str] = []
        for _, row in sub.iterrows():
            if current_start is None or row["qtl_start"] > current_end:
                if current_start is not None:
                    hotspot_rows.append(
                        {
                            "hotspot_id": f"HS{hotspot_id:04d}",
                            "chrom": chrom,
                            "start": current_start,
                            "end": current_end,
                            "trait_count": len(sorted(set(current_traits))),
                            "trait_ids": ",".join(sorted(set(current_traits))),
                        }
                    )
                    hotspot_id += 1
                current_start = int(row["qtl_start"])
                current_end = int(row["qtl_end"])
                current_traits = extract_traits(row)
            else:
                current_end = max(current_end, int(row["qtl_end"]))
                current_traits.extend(extract_traits(row))
        if current_start is not None:
            hotspot_rows.append(
                {
                    "hotspot_id": f"HS{hotspot_id:04d}",
                    "chrom": chrom,
                    "start": current_start,
                    "end": current_end,
                    "trait_count": len(sorted(set(current_traits))),
                    "trait_ids": ",".join(sorted(set(current_traits))),
                }
            )
            hotspot_id += 1
    return pd.DataFrame(hotspot_rows)


def add_representative_genes(candidate: pd.DataFrame, qtl_regions: pd.DataFrame) -> pd.DataFrame:
    if candidate.empty:
        return pd.DataFrame(columns=["qtl_id", "representative_gene_id", "distance_to_lead", "reason"])
    qtl_lead = qtl_regions.set_index("qtl_id")["lead_pos"].to_dict()
    rep = candidate.copy()
    rep["lead_pos"] = rep["qtl_id"].map(qtl_lead)
    rep["distance_to_lead"] = rep.apply(
        lambda row: 0
        if row["start"] <= row["lead_pos"] <= row["end"]
        else min(abs(row["start"] - row["lead_pos"]), abs(row["end"] - row["lead_pos"])),
        axis=1,
    )
    rep["inside_lead_gene"] = rep.apply(lambda row: row["start"] <= row["lead_pos"] <= row["end"], axis=1)
    rep = rep.sort_values(["qtl_id", "inside_lead_gene", "distance_to_lead", "gene_id"], ascending=[True, False, True, True])
    rep = rep.drop_duplicates(subset=["qtl_id"]).copy()
    rep["reason"] = rep["inside_lead_gene"].map(lambda x: "lead_snp_inside_gene" if x else "nearest_gene_to_lead_snp")
    rep = rep.rename(columns={"gene_id": "representative_gene_id"})
    return rep


def assign_significant_to_qtl(sig: pd.DataFrame, qtl_loci: pd.DataFrame) -> pd.DataFrame:
    assigned_parts = []
    for chrom, sig_sub in sig.groupby("chrom"):
        loci = qtl_loci[qtl_loci["chrom"] == chrom].sort_values("qtl_start").reset_index(drop=True)
        if loci.empty:
            continue
        starts = loci["qtl_start"].to_numpy()
        ends = loci["qtl_end"].to_numpy()
        qtl_ids = loci["qtl_id"].to_numpy()
        lead_pos = loci["lead_pos"].to_numpy()
        positions = sig_sub["pos"].to_numpy()
        idx = np.searchsorted(starts, positions, side="right") - 1
        valid = idx >= 0
        idx_valid = idx[valid]
        pos_valid = positions[valid]
        within = pos_valid <= ends[idx_valid]
        if not within.any():
            continue
        sub = sig_sub.iloc[np.where(valid)[0][within]].copy()
        chosen = idx_valid[within]
        sub["qtl_id"] = qtl_ids[chosen]
        sub["distance_to_lead"] = np.abs(sub["pos"].to_numpy() - lead_pos[chosen])
        assigned_parts.append(sub)
    if not assigned_parts:
        return pd.DataFrame(columns=list(sig.columns) + ["qtl_id", "distance_to_lead"])
    return pd.concat(assigned_parts, ignore_index=True)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lead = pd.read_csv(args.lead_snps, sep="\t")
    rename_map = {}
    if "trait_count" in lead.columns:
        rename_map["trait_count"] = "lead_trait_count"
    if "trait_ids" in lead.columns:
        rename_map["trait_ids"] = "lead_trait_ids"
    if rename_map:
        lead = lead.rename(columns=rename_map)
    summaries = pd.read_csv(args.trait_summaries, sep="\t") if args.trait_summaries and Path(args.trait_summaries).exists() else pd.DataFrame(columns=["trait_id", "trait_level"])
    annotation_available = bool(args.annotation_dir and Path(args.annotation_dir).exists())
    if annotation_available:
        annotation_dir = Path(args.annotation_dir)
        gene_meta = pd.read_csv(annotation_dir / "gene_metadata.tsv", sep="\t")
        if args.gene_annotation_file and Path(args.gene_annotation_file).exists():
            extra_annot = pd.read_csv(args.gene_annotation_file, sep="\t")
            if "gene_id" not in extra_annot.columns:
                raise ValueError("Gene annotation file must contain a gene_id column.")
            gene_meta = gene_meta.merge(extra_annot, on="gene_id", how="left", suffixes=("", "_extra"))
        genes = gene_meta.loc[:, [col for col in gene_meta.columns if col in {"chrom", "start", "end", "gene_id", "nr_annotation", "pfam"} or col not in {"bed_start", "promoter_start", "promoter_end"}]].copy()
        genes["chrom_norm"] = genes["chrom"].map(normalize_chrom_name)
        promoters = pd.read_csv(annotation_dir / "promoters_2kb.bed", sep="\t", header=None, names=["chrom", "start0", "end", "gene_id", "strand"])
        promoters["start"] = promoters["start0"] + 1
        promoters["chrom_norm"] = promoters["chrom"].map(normalize_chrom_name)
        downstream = pd.read_csv(annotation_dir / "downstream_2kb.bed", sep="\t", header=None, names=["chrom", "start0", "end", "gene_id", "strand"])
        downstream["start"] = downstream["start0"] + 1
        downstream["chrom_norm"] = downstream["chrom"].map(normalize_chrom_name)
        exons = pd.read_csv(annotation_dir / "exons.bed", sep="\t", header=None, names=["chrom", "start0", "end", "gene_id", "strand"])
        exons["start"] = exons["start0"] + 1
        exons["chrom_norm"] = exons["chrom"].map(normalize_chrom_name)
        cds = pd.read_csv(annotation_dir / "cds.bed", sep="\t", header=None, names=["chrom", "start0", "end", "gene_id", "strand"])
        cds["start"] = cds["start0"] + 1
        cds["chrom_norm"] = cds["chrom"].map(normalize_chrom_name)
    else:
        gene_meta = pd.DataFrame(columns=["gene_id"])
        genes = pd.DataFrame(columns=["chrom", "start", "end", "gene_id", "nr_annotation", "pfam", "chrom_norm"])
        promoters = pd.DataFrame(columns=["chrom", "start", "end", "gene_id", "chrom_norm"])
        downstream = pd.DataFrame(columns=["chrom", "start", "end", "gene_id", "chrom_norm"])
        exons = pd.DataFrame(columns=["chrom", "start", "end", "gene_id", "chrom_norm"])
        cds = pd.DataFrame(columns=["chrom", "start", "end", "gene_id", "chrom_norm"])

    ld_dir = Path(args.ld_dir)
    ld_records = []
    for _, row in lead.iterrows():
        ld_file = ld_dir / f"{str(row['snp_id']).replace(':', '_')}.ld"
        ld_records.append(load_ld_region(ld_file, row))
    ld_regions = pd.DataFrame(ld_records)
    qtl_loci = lead.merge(ld_regions, left_on="snp_id", right_on="lead_snp", how="left")
    if "chrom_x" in qtl_loci.columns:
        qtl_loci = qtl_loci.rename(columns={"chrom_x": "chrom"})
    if "chrom_y" in qtl_loci.columns:
        qtl_loci = qtl_loci.drop(columns=["chrom_y"])
    qtl_loci["qtl_id"] = [f"QTL{idx+1:05d}" for idx in range(qtl_loci.shape[0])]
    qtl_loci["qtl_size_kb"] = qtl_loci["qtl_size_bp"] / 1000
    qtl_loci["chrom_norm"] = qtl_loci["chrom"].map(normalize_chrom_name)
    qtl_regions = qtl_loci.copy()
    if "lead_trait_count" in qtl_regions.columns:
        qtl_regions["trait_count"] = qtl_regions["lead_trait_count"].fillna(0).astype(int)
    else:
        qtl_regions["trait_count"] = 0
    if "lead_trait_ids" in qtl_regions.columns:
        qtl_regions["trait_ids"] = qtl_regions["lead_trait_ids"].fillna("")
    else:
        qtl_regions["trait_ids"] = ""
    qtl_regions["lipid_count"] = qtl_regions["trait_ids"].map(
        lambda text: sum(item.startswith("molecule__") for item in str(text).split(",") if item)
    )
    qtl_regions.to_csv(outdir / "qtl_regions.tsv", sep="\t", index=False)

    trait_level_map = summaries.set_index("trait_id")["trait_level"].to_dict()
    qtl_trait_rows = []
    for _, row in qtl_regions.iterrows():
        for trait_id in [item for item in str(row.get("trait_ids", "")).split(",") if item]:
            qtl_trait_rows.append(
                {
                    "qtl_id": row["qtl_id"],
                    "trait_id": trait_id,
                    "trait_level": trait_level_map.get(trait_id, "unknown"),
                    "min_p": row.get("min_p", np.nan),
                    "significant_snp_count": 1,
                }
            )
    qtl_trait = pd.DataFrame(qtl_trait_rows, columns=["qtl_id", "trait_id", "trait_level", "min_p", "significant_snp_count"])
    qtl_trait.to_csv(outdir / "qtl_trait_membership.tsv", sep="\t", index=False)

    if not qtl_regions.empty:
        hotspot = build_hotspots(qtl_regions[["chrom", "qtl_start", "qtl_end", "trait_ids"]].copy())
        hotspot_cutoff = max(5, int(np.ceil(qtl_regions["trait_count"].quantile(0.95))))
        hotspot["hotspot_definition"] = f"merged overlapping QTL intervals; reference trait-count cutoff among loci = {hotspot_cutoff}"
    else:
        hotspot = pd.DataFrame(columns=["hotspot_id", "chrom", "start", "end", "trait_count", "trait_ids", "hotspot_definition"])
    hotspot.to_csv(outdir / "qtl_hotspots.tsv", sep="\t", index=False)

    candidate_rows = []
    if annotation_available:
        for _, row in qtl_regions.iterrows():
            chrom_norm = row["chrom_norm"]
            overlap = genes[(genes["chrom_norm"] == chrom_norm) & (genes["end"] >= row["qtl_start"]) & (genes["start"] <= row["qtl_end"])].copy()
            if overlap.empty:
                continue
            overlap["qtl_id"] = row["qtl_id"]
            overlap["lead_snp"] = row["snp_id"]
            overlap["trait_ids"] = row.get("trait_ids", "")
            overlap["trait_count"] = row.get("trait_count", 0)
            overlap["lead_pos"] = row["lead_pos"]
            candidate_rows.append(overlap)
    candidate = pd.concat(candidate_rows, ignore_index=True) if candidate_rows else pd.DataFrame(columns=["chrom", "start", "end", "gene_id", "nr_annotation", "pfam", "qtl_id", "lead_snp", "trait_ids", "trait_count", "lead_pos"])
    candidate.to_csv(outdir / "candidate_genes.tsv", sep="\t", index=False)
    representative = add_representative_genes(candidate, qtl_regions)
    representative.to_csv(outdir / "qtl_representative_genes.tsv", sep="\t", index=False)

    if annotation_available:
        unique_sig_path = outdir / "unique_significant_snps.tsv"
        if unique_sig_path.exists():
            unique_sig = pd.read_csv(unique_sig_path, sep="\t").loc[:, ["snp_id", "chrom", "pos"]].drop_duplicates().copy()
        else:
            unique_sig = pd.read_csv(args.master_significant, sep="\t", usecols=["snp_id", "chrom", "pos"]).drop_duplicates().copy()
        unique_sig["trait_id"] = "multi_trait"
        snp_anno = annotate_snps(unique_sig, genes, promoters[["chrom", "start", "end", "gene_id"]], downstream[["chrom", "start", "end", "gene_id"]], exons[["chrom", "start", "end", "gene_id"]], cds[["chrom", "start", "end", "gene_id"]])
        snp_anno.to_csv(outdir / "significant_snp_annotation.tsv.gz", sep="\t", index=False)
    else:
        pd.DataFrame(columns=["snp_id", "trait_id", "chrom", "pos", "region_type", "gene_id"]).to_csv(outdir / "significant_snp_annotation.tsv.gz", sep="\t", index=False)

    qtl_trait_counts = qtl_regions.loc[:, ["qtl_id", "trait_count", "trait_ids", "lipid_count"]].copy()
    qtl_trait_counts.to_csv(outdir / "qtl_trait_counts.tsv", sep="\t", index=False)
    qtl_size_stats = qtl_regions["qtl_size_kb"].describe().to_frame().reset_index().rename(columns={"index": "stat", "qtl_size_kb": "value"})
    qtl_size_stats.to_csv(outdir / "qtl_size_stats.tsv", sep="\t", index=False)

    if args.module_membership and Path(args.module_membership).exists():
        modules = pd.read_csv(args.module_membership, sep="\t")
        module_qtl = qtl_trait.merge(modules[["trait_id", "module_id"]], on="trait_id", how="left")
        module_qtl = module_qtl.groupby(["module_id", "qtl_id"]).agg(trait_count=("trait_id", "nunique")).reset_index()
        module_qtl.to_csv(outdir / "module_qtl_summary.tsv", sep="\t", index=False)

    (outdir / "non_synonymous_status.txt").write_text(
        "Non-synonymous SNP identification was not performed because the supplied genotype files encode alleles as 1/2 rather than A/C/G/T bases, which prevents codon-level consequence inference.\n"
    )
    (outdir / "pathway_note.txt").write_text(
        "Pathway enrichment here is keyword-driven using nr_annotation/pfam because no direct gene-to-KEGG mapping file was provided in the workspace.\n"
    )

    check = {
        "qtl_count": int(qtl_regions.shape[0]),
        "hotspot_count": int(hotspot.shape[0]),
        "candidate_gene_count": int(candidate["gene_id"].nunique()) if not candidate.empty else 0,
        "note": "QTL regions are defined from PLINK LD windows around lead SNPs (R2 >= 0.6, 500 kb search window).",
    }
    (outdir / "qtl_check.json").write_text(json.dumps(check, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
