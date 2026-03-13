#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


CHAIN_RE = re.compile(r"(\d+):(\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare lipid phenotype matrices for GWAS.")
    parser.add_argument("--id-file", required=True)
    parser.add_argument("--lipid-file", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def inverse_normal_transform(values: pd.Series) -> pd.Series:
    observed = values.dropna()
    if observed.empty:
        return pd.Series(np.nan, index=values.index)
    if observed.nunique() <= 1:
        return pd.Series(0.0, index=values.index)
    ranks = rankdata(observed.values, method="average")
    transformed = norm.ppf((ranks - 0.5) / len(observed))
    result = pd.Series(np.nan, index=values.index, dtype=float)
    result.loc[observed.index] = transformed
    if result.isna().any():
        result = result.fillna(result.mean())
    return result


def classify_superclass(lipid_class: str) -> str:
    glycerolipids = {"TG", "DG", "MG", "DGDG", "MGDG", "DGMG", "MGMG", "SQDG"}
    glycerophospholipids = {"PC", "PE", "PG", "PI", "PS", "PA", "CL", "LPC", "LPE", "LPG", "LPI", "LPS"}
    sphingolipids = {
        "Cer", "CerP", "SPH", "Hex1Cer", "Hex2Cer", "Hex3Cer", "Hex1SPH",
        "CerG2GNAc1", "CerG3GNAc1", "GM2", "GM3", "GD1a", "GD2", "GD3",
    }
    sterol_related = {"ChE", "Co", "WE"}
    fatty_acyls = {"AcCa", "FAHFA", "FA", "CAR"}

    if lipid_class in glycerolipids:
        return "Glycerolipid"
    if lipid_class in glycerophospholipids:
        return "Glycerophospholipid"
    if lipid_class in sphingolipids:
        return "Sphingolipid"
    if lipid_class in sterol_related:
        return "Sterol_Or_Wax"
    if lipid_class in fatty_acyls:
        return "Fatty_Acyl"
    return "Other"


def classify_role(lipid_class: str) -> str:
    storage = {"TG", "DG", "MG", "WE", "ChE"}
    membrane = {"PC", "PE", "PG", "PI", "PS", "PA", "CL", "LPC", "LPE", "LPG", "LPI", "LPS", "DGDG", "MGDG", "SQDG", "MGMG", "DGMG"}
    signaling = {"CerP", "SPH", "Co", "AcCa"}
    if lipid_class in storage:
        return "Storage"
    if lipid_class in membrane:
        return "Membrane"
    if lipid_class in signaling:
        return "Signaling"
    if lipid_class.startswith("Cer") or "Cer" in lipid_class or lipid_class.startswith(("GM", "GD", "Hex")):
        return "Sphingo_Glyco"
    return "Other"


def parse_chain_stats(row: pd.Series) -> tuple[float | None, float | None, int]:
    values: list[tuple[int, int]] = []
    for col in ["FA1", "FA2", "FA3", "FA4"]:
        raw = str(row.get(col, "")).strip()
        if raw and raw != "N/A":
            values.extend((int(a), int(b)) for a, b in CHAIN_RE.findall(raw))
    if not values:
        raw = str(row.get("FattyAcid", "")).strip()
        values.extend((int(a), int(b)) for a, b in CHAIN_RE.findall(raw))
    if not values:
        raw = str(row.get("LipidIon", "")).strip()
        values.extend((int(a), int(b)) for a, b in CHAIN_RE.findall(raw))
    if not values:
        return None, None, 0
    total_c = sum(a for a, _ in values)
    total_db = sum(b for _, b in values)
    return float(total_c), float(total_db), len(values)


def chain_length_bin(total_c: float | None) -> str:
    if total_c is None:
        return "Unknown"
    if total_c <= 34:
        return "C_le_34"
    if total_c <= 40:
        return "C35_40"
    if total_c <= 46:
        return "C41_46"
    return "C_ge_47"


def unsaturation_bin(total_db: float | None) -> str:
    if total_db is None:
        return "Unknown"
    if total_db == 0:
        return "DB0"
    if total_db <= 2:
        return "DB1_2"
    if total_db <= 4:
        return "DB3_4"
    return "DB_ge_5"


def build_group_traits(
    row_means: pd.DataFrame,
    metadata: pd.DataFrame,
    group_col: str,
    prefix: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records = []
    trait_meta = []
    for group_name, sub_meta in metadata.groupby(group_col, dropna=False):
        if pd.isna(group_name):
            continue
        row_ids = sub_meta["trait_source_id"].tolist()
        trait_id = f"{prefix}__{group_name}"
        series = row_means.loc[row_ids].mean(axis=0, skipna=True)
        series.name = trait_id
        records.append(series)
        trait_meta.append(
            {
                "trait_id": trait_id,
                "trait_level": prefix,
                "display_name": group_name,
                "member_count": len(row_ids),
                "source_group": group_name,
            }
        )
    if not records:
        return pd.DataFrame(), pd.DataFrame()
    out = pd.DataFrame(records)
    out.index = [x["trait_id"] for x in trait_meta]
    return out, pd.DataFrame(trait_meta)


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "by_level").mkdir(exist_ok=True)
    (outdir / "checks").mkdir(exist_ok=True)

    ids = [line.strip() for line in Path(args.id_file).read_text().splitlines() if line.strip()]
    lipid = pd.read_csv(args.lipid_file, sep="\t", dtype=str)
    sample_cols = lipid.columns[17:].tolist()

    group_to_cols: dict[str, list[str]] = defaultdict(list)
    for col in sample_cols:
        if "_" not in col:
            continue
        prefix, suffix = col.rsplit("_", 1)
        if suffix not in {"1", "2", "3"}:
            continue
        group_to_cols[prefix].append(col)

    errors: list[str] = []
    for sample_id in ids:
        cols = group_to_cols.get(sample_id, [])
        if len(cols) != 3:
            errors.append(f"{sample_id} has {len(cols)} replicate columns")
    if errors:
        raise ValueError("Replicate validation failed: " + "; ".join(errors[:10]))

    numeric_values = lipid[sample_cols].apply(pd.to_numeric, errors="coerce")
    row_means = pd.DataFrame(
        {
            sample_id: numeric_values[group_to_cols[sample_id]].mean(axis=1, skipna=True).to_numpy()
            for sample_id in ids
        },
        index=lipid.index.astype(str),
    )

    metadata = lipid[[
        "name", "DIFF", "LipidIon", "LipidGroup", "Class", "FattyAcid", "FA1", "FA2", "FA3", "FA4",
        "IonFormula", "CalMz", "RT_(min)", "QCRSD", "Fold Change", "P_value", "VIP"
    ]].copy()
    metadata["trait_source_id"] = lipid.index.astype(str)
    metadata["trait_id"] = metadata["name"].astype(str).map(lambda x: f"molecule__{x}")
    metadata["display_name"] = metadata["LipidIon"].fillna(metadata["name"])
    metadata["trait_level"] = "molecule"

    chain_stats = metadata.apply(parse_chain_stats, axis=1, result_type="expand")
    chain_stats.columns = ["total_carbon", "total_double_bond", "chain_count"]
    metadata = pd.concat([metadata, chain_stats], axis=1)
    metadata["superclass"] = metadata["Class"].astype(str).map(classify_superclass)
    metadata["biological_role"] = metadata["Class"].astype(str).map(classify_role)
    metadata["chain_length_bin"] = metadata["total_carbon"].map(chain_length_bin)
    metadata["unsaturation_bin"] = metadata["total_double_bond"].map(unsaturation_bin)
    metadata["chain_count_bin"] = metadata["chain_count"].map(
        lambda x: "Unknown" if pd.isna(x) or x == 0 else f"{int(x)}_chains"
    )

    row_means.index = metadata["trait_source_id"]

    molecule_matrix = row_means.copy()
    molecule_matrix.index = metadata["trait_id"].tolist()

    level_tables: dict[str, pd.DataFrame] = {"molecule": molecule_matrix}
    level_meta = [metadata[[
        "trait_id", "trait_level", "display_name", "Class", "LipidGroup", "FattyAcid", "QCRSD",
        "Fold Change", "P_value", "VIP", "superclass", "biological_role", "chain_length_bin",
        "unsaturation_bin", "chain_count_bin", "total_carbon", "total_double_bond", "chain_count"
    ]].copy()]

    for group_col, prefix in [
        ("Class", "class"),
        ("superclass", "superclass"),
        ("biological_role", "role"),
        ("chain_length_bin", "chainlen"),
        ("unsaturation_bin", "unsat"),
        ("chain_count_bin", "chaincount"),
    ]:
        table, meta = build_group_traits(row_means, metadata, group_col, prefix)
        if not table.empty:
            level_tables[prefix] = table
            level_meta.append(meta)

    total_series = row_means.mean(axis=0, skipna=True)
    total_series.name = "total__all_lipids"
    total_table = pd.DataFrame([total_series])
    level_tables["total"] = total_table
    level_meta.append(
        pd.DataFrame(
            [
                {
                    "trait_id": "total__all_lipids",
                    "trait_level": "total",
                    "display_name": "All lipids mean abundance",
                    "member_count": row_means.shape[0],
                    "source_group": "all_lipids",
                }
            ]
        )
    )

    combined = pd.concat(level_tables.values(), axis=0)
    combined = combined.loc[~combined.index.duplicated()].copy()

    combined = combined[ids]
    combined.columns.name = "sample_id"

    raw_matrix = combined.T
    raw_matrix.index.name = "IID"
    raw_matrix.insert(0, "IID", raw_matrix.index)
    raw_matrix.insert(0, "FID", raw_matrix.index)

    int_only = combined.T.apply(inverse_normal_transform, axis=0)
    int_only.index.name = "IID"
    int_matrix = int_only.copy()
    int_matrix.insert(0, "IID", int_matrix.index)
    int_matrix.insert(0, "FID", int_matrix.index)

    missing_rate = raw_matrix.iloc[:, 2:].isna().mean(axis=0)
    raw_matrix.iloc[:, 2:] = raw_matrix.iloc[:, 2:].fillna(raw_matrix.iloc[:, 2:].median(axis=0))
    int_matrix.iloc[:, 2:] = int_matrix.iloc[:, 2:].fillna(0.0)

    trait_meta = pd.concat(level_meta, ignore_index=True, sort=False)
    trait_meta["missing_rate"] = trait_meta["trait_id"].map(missing_rate.to_dict()).fillna(0.0)
    trait_meta["sample_count"] = len(ids)
    trait_meta = trait_meta.drop_duplicates(subset=["trait_id"]).sort_values(["trait_level", "trait_id"])

    raw_matrix.to_csv(outdir / "all_traits_raw.tsv", sep="\t", index=False)
    int_matrix.to_csv(outdir / "all_traits_int.tsv", sep="\t", index=False)
    trait_meta.to_csv(outdir / "trait_metadata.tsv", sep="\t", index=False)

    summary_records = []
    for level_name, table in level_tables.items():
        raw_subset = raw_matrix[["FID", "IID"] + table.index.tolist()]
        int_subset = int_matrix[["FID", "IID"] + table.index.tolist()]
        raw_subset.to_csv(outdir / "by_level" / f"{level_name}_raw.tsv", sep="\t", index=False)
        int_subset.to_csv(outdir / "by_level" / f"{level_name}_int.tsv", sep="\t", index=False)
        summary_records.append(
            {
                "trait_level": level_name,
                "trait_count": table.shape[0],
                "sample_count": table.shape[1],
            }
        )

    pd.DataFrame(summary_records).to_csv(outdir / "checks" / "trait_level_counts.tsv", sep="\t", index=False)

    check_note = {
        "ids_count": len(ids),
        "sample_columns": len(sample_cols),
        "replicate_rule": "Each sample_id must have exactly 3 replicate columns ending with _1/_2/_3.",
        "notes": [
            "Replicate means are computed by sample_id prefix, not by assuming adjacent columns.",
            "Missing phenotype values are median-imputed for raw matrices and zero-imputed after INT.",
            "Additional phenotype levels include class, superclass, biological role, chain length bin, unsaturation bin, chain-count bin, and total lipid mean.",
        ],
    }
    (outdir / "checks" / "prepare_phenotypes_check.json").write_text(
        json.dumps(check_note, indent=2, ensure_ascii=False) + "\n"
    )


if __name__ == "__main__":
    main()
