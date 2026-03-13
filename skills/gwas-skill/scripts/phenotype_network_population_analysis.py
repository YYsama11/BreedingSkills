#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.vq import kmeans2
from scipy.stats import f_oneway, pearsonr


plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.dpi"] = 150
sns.set_theme(style="whitegrid")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Phenotype, network, and population analyses for lipid GWAS.")
    parser.add_argument("--raw-phenotypes", required=True)
    parser.add_argument("--int-phenotypes", required=True)
    parser.add_argument("--trait-meta", required=True)
    parser.add_argument("--trait-summaries", required=True)
    parser.add_argument("--pca", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def run_pca(matrix: pd.DataFrame, n_components: int = 5) -> tuple[np.ndarray, np.ndarray]:
    values = matrix.to_numpy(dtype=float)
    values = values - values.mean(axis=0, keepdims=True)
    values = values / np.where(values.std(axis=0, keepdims=True) == 0, 1, values.std(axis=0, keepdims=True))
    u, s, _ = np.linalg.svd(values, full_matrices=False)
    scores = u[:, :n_components] * s[:n_components]
    explained = (s ** 2) / np.sum(s ** 2)
    return scores, explained[:n_components]


def main() -> None:
    args = parse_args()
    outdir = Path(args.outdir)
    fig_dir = outdir / "figures"
    for directory in [outdir, fig_dir]:
        directory.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.raw_phenotypes, sep="\t")
    int_df = pd.read_csv(args.int_phenotypes, sep="\t")
    trait_meta = pd.read_csv(args.trait_meta, sep="\t")
    summaries = pd.read_csv(args.trait_summaries, sep="\t")
    pca = pd.read_csv(args.pca, sep=r"\s+", header=None)
    cov = pd.read_csv(args.covariates, sep="\t")

    raw = raw.set_index("IID")
    int_df = int_df.set_index("IID")
    molecule_traits = trait_meta.loc[trait_meta["trait_level"] == "molecule", "trait_id"].tolist()
    aggregate_traits = trait_meta.loc[trait_meta["trait_level"].isin(["class", "superclass", "role", "total"]), "trait_id"].tolist()

    pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pca.shape[1] - 1)]
    pca = pca.set_index("IID")
    subgroup_input = pca[["PC1", "PC2"]].to_numpy()
    _, labels = kmeans2(subgroup_input, 3, minit="points")
    pca["subgroup"] = [f"Group{label+1}" for label in labels]
    pca.to_csv(outdir / "sample_groups.tsv", sep="\t")
    for subgroup, sub_df in pca.groupby("subgroup"):
        keep = sub_df.reset_index()[["IID", "IID"]]
        keep.to_csv(outdir / f"{subgroup}.keep.txt", sep="\t", index=False, header=False)

    total_trait = "total__all_lipids"
    total_values = raw[total_trait]
    fig, ax = plt.subplots(figsize=(6, 5))
    palette = sns.color_palette("RdYlGn", n_colors=pca["subgroup"].nunique())
    sns.scatterplot(data=pca, x="PC1", y="PC2", hue="subgroup", size=total_values.loc[pca.index], palette=palette, ax=ax)
    ax.set_title("Population structure vs total lipid abundance")
    fig.tight_layout()
    fig.savefig(fig_dir / "population_structure_total_lipid.png")
    plt.close(fig)

    # Phenotype PCA
    agg_scores, agg_var = run_pca(int_df[aggregate_traits], n_components=3)
    agg_pca = pd.DataFrame(agg_scores, index=int_df.index, columns=["PC1", "PC2", "PC3"]).join(pca[["subgroup"]], how="left")
    agg_pca.to_csv(outdir / "aggregate_trait_pca_scores.tsv", sep="\t")
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.scatterplot(data=agg_pca, x="PC1", y="PC2", hue="subgroup", palette=palette, ax=ax)
    ax.set_title(f"Phenotype PCA (aggregate traits)\nVar={agg_var[0]:.2%},{agg_var[1]:.2%}")
    fig.tight_layout()
    fig.savefig(fig_dir / "aggregate_trait_pca.png")
    plt.close(fig)

    # Clustering heatmap on aggregate traits
    heatmap_df = raw[aggregate_traits].copy()
    heatmap_df.index.name = "sample_id"
    clustermap = sns.clustermap(heatmap_df, cmap="RdYlGn", figsize=(10, 10), col_cluster=True, row_cluster=True)
    clustermap.fig.savefig(fig_dir / "aggregate_trait_clustermap.png")
    plt.close(clustermap.fig)

    # Group differences
    diff_records = []
    for trait in aggregate_traits:
        groups = []
        for subgroup, sub in pca.groupby("subgroup"):
            groups.append(raw.loc[sub.index, trait].astype(float).to_numpy())
        stat = f_oneway(*groups)
        means = raw.loc[pca.index, [trait]].join(pca[["subgroup"]]).groupby("subgroup")[trait].mean().to_dict()
        diff_records.append({"trait_id": trait, "p_value": stat.pvalue, **{f"mean_{k}": v for k, v in means.items()}})
    diff_df = pd.DataFrame(diff_records).sort_values("p_value")
    diff_df.to_csv(outdir / "subgroup_trait_differences.tsv", sep="\t", index=False)

    top_diff = diff_df.head(20)["trait_id"].tolist()
    if top_diff:
        plot_df = raw.loc[pca.index, top_diff].join(pca[["subgroup"]])
        melted = plot_df.reset_index().melt(id_vars=["IID", "subgroup"], var_name="trait_id", value_name="value")
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.boxplot(data=melted, x="trait_id", y="value", hue="subgroup", palette=palette, ax=ax)
        ax.tick_params(axis="x", rotation=90)
        ax.set_title("Top subgroup-differentiated aggregate lipid traits")
        fig.tight_layout()
        fig.savefig(fig_dir / "subgroup_trait_differences_boxplot.png")
        plt.close(fig)

    # PC-lipid correlation
    corr_records = []
    for pc in ["PC1", "PC2", "PC3", "PC4", "PC5"]:
        if pc not in pca.columns:
            continue
        for trait in aggregate_traits:
            r, p_value = pearsonr(pca.loc[int_df.index, pc], int_df[trait])
            corr_records.append({"pc": pc, "trait_id": trait, "correlation": r, "p_value": p_value})
    corr_df = pd.DataFrame(corr_records).sort_values("p_value")
    corr_df.to_csv(outdir / "population_structure_lipid_correlations.tsv", sep="\t", index=False)

    # Network analysis
    molecule_matrix = int_df[molecule_traits]
    corr = molecule_matrix.corr(method="pearson")
    corr.to_csv(outdir / "molecule_correlation_matrix.tsv.gz", sep="\t")
    edges = []
    cols = corr.columns.tolist()
    for i, col_i in enumerate(cols):
        values = corr.iloc[i, i + 1 :].to_numpy()
        idx = np.where(np.abs(values) >= 0.6)[0]
        for offset in idx:
            j = i + 1 + offset
            edges.append((col_i, cols[j], float(values[offset])))
    edge_df = pd.DataFrame(edges, columns=["source", "target", "correlation"])
    edge_df.to_csv(outdir / "lipid_network_edges.tsv.gz", sep="\t", index=False)

    graph = nx.Graph()
    for source, target, value in edges:
        graph.add_edge(source, target, weight=value)
    graph.add_nodes_from(molecule_traits)

    degree = dict(graph.degree())
    degree_centrality = nx.degree_centrality(graph)
    communities = list(nx.algorithms.community.greedy_modularity_communities(graph))
    module_map = {}
    for idx, community in enumerate(communities, start=1):
        for node in community:
            module_map[node] = f"M{idx:03d}"

    nodes = pd.DataFrame({"trait_id": list(graph.nodes())})
    nodes["degree"] = nodes["trait_id"].map(degree).fillna(0)
    nodes["degree_centrality"] = nodes["trait_id"].map(degree_centrality).fillna(0)
    nodes["module_id"] = nodes["trait_id"].map(module_map).fillna("M000")
    nodes = nodes.merge(trait_meta[["trait_id", "Class", "display_name"]], on="trait_id", how="left")
    nodes = nodes.merge(summaries[["trait_id", "significant_snp_count", "top_p"]], on="trait_id", how="left")
    nodes["hub_flag"] = nodes["degree"] >= nodes["degree"].quantile(0.95)
    nodes.to_csv(outdir / "lipid_network_nodes.tsv", sep="\t", index=False)
    nodes[["trait_id", "module_id"]].to_csv(outdir / "lipid_module_membership.tsv", sep="\t", index=False)

    top_nodes = nodes.sort_values(["degree", "significant_snp_count"], ascending=[False, False]).head(150)["trait_id"].tolist()
    subgraph = graph.subgraph(top_nodes).copy()
    if subgraph.number_of_nodes() > 0:
        pos = nx.spring_layout(subgraph, seed=42)
        fig, ax = plt.subplots(figsize=(10, 10))
        node_sizes = [100 + subgraph.degree(node) * 20 for node in subgraph.nodes()]
        node_colors = [nodes.set_index("trait_id").loc[node, "significant_snp_count"] if node in nodes["trait_id"].values else 0 for node in subgraph.nodes()]
        nx.draw_networkx_edges(subgraph, pos, alpha=0.15, ax=ax)
        nx.draw_networkx_nodes(subgraph, pos, node_size=node_sizes, node_color=node_colors, cmap="RdYlGn", ax=ax)
        ax.set_title("Hub lipid network with GWAS signal overlay")
        ax.axis("off")
        fig.tight_layout()
        fig.savefig(fig_dir / "hub_lipid_network_gwas_overlay.png")
        plt.close(fig)

    h2_summary = summaries[["trait_id", "trait_level", "pseudo_heritability"]].copy()
    h2_summary.to_csv(outdir / "pseudo_heritability_summary.tsv", sep="\t", index=False)

    check = {
        "subgroup_count": int(pca["subgroup"].nunique()),
        "network_node_count": int(nodes.shape[0]),
        "network_edge_count": int(edge_df.shape[0]),
        "module_count": int(len(set(module_map.values()))),
    }
    (outdir / "phenotype_network_check.json").write_text(json.dumps(check, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
