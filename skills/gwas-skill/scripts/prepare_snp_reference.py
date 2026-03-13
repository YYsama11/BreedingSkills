#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare SNP reference arrays for GWAS plotting and summaries.")
    parser.add_argument("--bim", required=True)
    parser.add_argument("--fai", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fai = pd.read_csv(args.fai, sep="\t", header=None, names=["chrom", "length", "offset", "linebases", "linewidth"])
    chrom_order = [f"Chr{int(str(chrom).replace('Chr', ''))}" for chrom in fai["chrom"]]
    chrom_lengths = {int(str(chrom).replace("Chr", "")): int(length) for chrom, length in zip(fai["chrom"], fai["length"])}

    chrom_offsets = {}
    cumulative = 0
    tick_records = []
    for chrom in sorted(chrom_lengths):
        chrom_offsets[chrom] = cumulative
        tick_records.append({"chrom": chrom, "label": f"Chr{chrom}", "midpoint": cumulative + chrom_lengths[chrom] / 2})
        cumulative += chrom_lengths[chrom]

    bim = pd.read_csv(
        args.bim,
        sep="\t",
        header=None,
        names=["chrom", "snp_id", "cm", "pos", "a1", "a2"],
        dtype={"chrom": int, "snp_id": str, "cm": float, "pos": int, "a1": str, "a2": str},
    )
    cum_pos = bim["chrom"].map(chrom_offsets).to_numpy(dtype=np.int64) + bim["pos"].to_numpy(dtype=np.int64)

    total_length = cumulative
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
        for chrom in sorted(chrom_offsets):
            start = chrom_offsets[chrom]
            end = start + chrom_lengths[chrom]
            if start <= mid <= end:
                found = chrom
                break
        bin_chrom.append(found if found is not None else max(chrom_offsets))
    bin_table["chrom"] = bin_chrom
    bin_table["color_group"] = bin_table["chrom"] % 2

    np.save(outdir / "chrom.npy", bim["chrom"].to_numpy(dtype=np.int16))
    np.save(outdir / "pos.npy", bim["pos"].to_numpy(dtype=np.int32))
    np.save(outdir / "cum_pos.npy", cum_pos)
    np.save(outdir / "bin_index.npy", bin_index)
    bim[["snp_id", "a1", "a2"]].to_csv(outdir / "snp_id_alleles.tsv.gz", sep="\t", index=False)
    pd.DataFrame(tick_records).to_csv(outdir / "chrom_ticks.tsv", sep="\t", index=False)
    bin_table.to_csv(outdir / "manhattan_bins.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"metric": "snp_count", "value": len(bim)}, {"metric": "bonferroni_threshold", "value": 0.05 / len(bim)}, {"metric": "suggestive_threshold", "value": 1.0 / len(bim)}]
    ).to_csv(outdir / "reference_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
