#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare SNP reference arrays for GWAS plotting and summaries.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--bim", required=False)
    parser.add_argument("--tped", required=False)
    parser.add_argument("--fai", required=False)
    args = parser.parse_args()
    if not args.bim and not args.tped:
        parser.error("One of --bim or --tped is required.")
    return args


def normalize_chrom_name(name: str) -> str:
    text = str(name).strip()
    lowered = text.lower()
    if lowered.startswith("chr"):
        lowered = lowered[3:]
    if lowered.isdigit():
        lowered = str(int(lowered))
    return lowered


def read_variant_table(args: argparse.Namespace) -> pd.DataFrame:
    if args.bim:
        return pd.read_csv(
            args.bim,
            sep=r"\s+",
            header=None,
            names=["chrom", "snp_id", "cm", "pos", "a1", "a2"],
            dtype={"chrom": str, "snp_id": str, "cm": float, "pos": int, "a1": str, "a2": str},
        )
    return pd.read_csv(
        args.tped,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["chrom", "snp_id", "cm", "pos"],
        dtype={"chrom": str, "snp_id": str, "cm": float, "pos": int},
    ).assign(a1="NA", a2="NA")


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = read_variant_table(args)
    if variants.empty:
        raise ValueError("No variants were read from the provided BIM/TPED input.")

    if args.fai:
        fai = pd.read_csv(args.fai, sep="\t", header=None, names=["chrom", "length", "offset", "linebases", "linewidth"], dtype={0: str})
        fai["chrom"] = fai["chrom"].astype(str)
        chrom_lengths = {chrom: int(length) for chrom, length in zip(fai["chrom"], fai["length"])}
        fai_norm = {}
        for chrom in fai["chrom"]:
            fai_norm.setdefault(normalize_chrom_name(chrom), []).append(chrom)
        chrom_map = {}
        unresolved = []
        for chrom in variants["chrom"].astype(str).unique():
            if chrom in chrom_lengths:
                chrom_map[chrom] = chrom
                continue
            norm = normalize_chrom_name(chrom)
            candidates = fai_norm.get(norm, [])
            if len(candidates) == 1:
                chrom_map[chrom] = candidates[0]
            else:
                unresolved.append(chrom)
        if unresolved:
            raise ValueError(f"Unable to map variant chromosome names to FAI chromosome names: {unresolved[:20]}")
        variants["chrom_name"] = variants["chrom"].map(chrom_map)
        chrom_order = fai["chrom"].tolist()
    else:
        variants["chrom_name"] = variants["chrom"].astype(str)
        chrom_order = list(dict.fromkeys(variants["chrom_name"].tolist()))
        chrom_lengths = variants.groupby("chrom_name")["pos"].max().astype(int).to_dict()

    chrom_offsets = {}
    cumulative = 0
    tick_records = []
    for chrom in chrom_order:
        chrom_offsets[chrom] = cumulative
        tick_records.append({"chrom": chrom, "label": chrom, "midpoint": cumulative + chrom_lengths[chrom] / 2})
        cumulative += chrom_lengths[chrom]

    cum_pos = variants["chrom_name"].map(chrom_offsets).to_numpy(dtype=np.int64) + variants["pos"].to_numpy(dtype=np.int64)
    total_length = max(cumulative, 1)
    bin_count = 20000
    bin_index = np.floor(cum_pos / total_length * bin_count).astype(np.int32)
    bin_index[bin_index == bin_count] = bin_count - 1

    bin_table = pd.DataFrame({"bin_index": np.arange(bin_count)})
    bin_table["bin_start"] = np.floor(bin_table["bin_index"] * total_length / bin_count).astype(int)
    bin_table["bin_end"] = np.floor((bin_table["bin_index"] + 1) * total_length / bin_count).astype(int)
    bin_table["bin_mid"] = ((bin_table["bin_start"] + bin_table["bin_end"]) / 2).astype(float)

    bin_chrom = []
    for mid in bin_table["bin_mid"]:
        found = None
        for chrom in chrom_order:
            start = chrom_offsets[chrom]
            end = start + chrom_lengths[chrom]
            if start <= mid <= end:
                found = chrom
                break
        bin_chrom.append(found if found is not None else chrom_order[-1])
    bin_table["chrom"] = bin_chrom
    chrom_to_group = {chrom: idx % 2 for idx, chrom in enumerate(chrom_order)}
    bin_table["color_group"] = bin_table["chrom"].map(chrom_to_group)

    max_len = max(len(chrom) for chrom in variants["chrom_name"].astype(str)) if not variants.empty else 1
    np.save(outdir / "chrom.npy", variants["chrom_name"].astype(str).to_numpy(dtype=f"<U{max_len}"))
    np.save(outdir / "pos.npy", variants["pos"].to_numpy(dtype=np.int32))
    np.save(outdir / "cum_pos.npy", cum_pos)
    np.save(outdir / "bin_index.npy", bin_index)
    variants[["snp_id", "a1", "a2"]].to_csv(outdir / "snp_id_alleles.tsv.gz", sep="\t", index=False)
    pd.DataFrame(tick_records).to_csv(outdir / "chrom_ticks.tsv", sep="\t", index=False)
    bin_table.to_csv(outdir / "manhattan_bins.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {"metric": "snp_count", "value": len(variants)},
            {"metric": "bonferroni_threshold", "value": 0.05 / len(variants)},
            {"metric": "suggestive_threshold", "value": 1.0 / len(variants)},
        ]
    ).to_csv(outdir / "reference_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
