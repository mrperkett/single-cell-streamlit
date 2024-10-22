import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import streamlit as st


def get_umi_count_gene_count_scatterplot(adata):
    fig, ax = plt.subplots()

    # sort a copy of adata.obs so that high percent mitochondrial barcodes are not buried
    df = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    df = df.sort_values(["pct_counts_mt"], ascending=True)

    # Create a scatterplot coloring by pct_counts_mt
    plt.scatter(
        x=df["total_counts"],
        y=df["n_genes_by_counts"],
        c=df["pct_counts_mt"],
        s=0.85,
        cmap="viridis",
    )
    plt.colorbar(label="% UMIs that are Mitochondrial")
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    ax.set_xlabel("UMIs / Barcode", fontsize=16)
    ax.set_ylabel("Genes Identified / Barcode", fontsize=16)
    ax.set_title("Mitochondrial QC", fontsize=18)

    fig = plt.gcf()
    return fig


def get_log_umi_count_log_gene_count_jointplot(adata):
    # inspired by: https://scanpy.readthedocs.io/en/1.10.x/generated/scanpy.pp.calculate_qc_metrics.html#scanpy.pp.calculate_qc_metrics
    fig, ax = plt.subplots()
    sns.jointplot(
        data=adata.obs,
        x="log1p_total_counts",
        y="log1p_n_genes_by_counts",
        kind="hex",
    )
    ax = plt.gca()
    ax.set_xlabel("log2(total_count + 1)", fontsize=16)
    ax.set_ylabel("log2(num_genes + 1)", fontsize=16)

    fig = plt.gcf()
    return fig


def get_umi_distribution_plot(adata):
    fig, ax = plt.subplots()
    total_umis_by_cell = np.array(adata.X.astype(int).sum(axis=1)).flatten()
    # sns.histplot(adata.obs["total_counts"], ax=ax)
    sns.histplot(total_umis_by_cell, ax=ax)
    ax.set_xlabel("UMI Count", fontsize=16)
    ax.set_ylabel("Barcode Count", fontsize=16)
    ax.set_title("Distribution for UMI Count / Barcode", fontsize=18)
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))

    fig = plt.gcf()
    return fig


def get_gene_distribution_plot(adata):
    fig, ax = plt.subplots()
    sns.histplot(adata.obs["n_genes_by_counts"], ax=ax)
    ax.set_xlabel("Gene Count", fontsize=16)
    ax.set_ylabel("Barcode Count", fontsize=16)
    ax.set_title("Distribution for Gene Count / Barcode", fontsize=18)

    fig = plt.gcf()
    return fig


def get_mitochondrial_umi_distribution_plot(adata):
    fig, ax = plt.subplots()
    sns.histplot(adata.obs["pct_counts_mt"], ax=ax)
    ax.set_xlabel("% UMIs that are Mitochondrial", fontsize=16)
    ax.set_ylabel("Barcode Count", fontsize=16)
    ax.set_title("Distribution of % Mitochondrial", fontsize=18)

    fig = plt.gcf()
    return fig


def display_umi_count_gene_count_scatterplot(adata):
    fig = get_umi_count_gene_count_scatterplot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")


def display_log_umi_count_log_gene_count_jointplot(adata):
    fig = get_log_umi_count_log_gene_count_jointplot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")


def display_get_umi_distribution_plot(adata):
    fig = get_umi_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")


def display_gene_distribution_plot(adata):
    fig = get_gene_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")


def display_mitochondrial_umi_distribution_plot(adata):
    fig = get_mitochondrial_umi_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def display_doublet_detection_info(adata, algorithm="Scrublet"):
    if algorithm == "Scrublet":
        display_doublet_detection_info_scrublet(adata)
    elif algorithm == "Vaeda":
        raise NotImplementedError("Vaeda algorithm not implemented")
    else:
        raise ValueError("algorithm '{algorithm}' not recgonized")


def display_doublet_detection_info_scrublet(adata):

    num_predicted_doublets = adata.obs["predicted_doublet"].sum()
    num_total = len(adata.obs)
    num_predicted_singlets = num_total - num_predicted_doublets

    data = [
        ["Singlets", num_predicted_singlets, num_predicted_singlets / num_total * 100.0],
        ["Doublets", num_predicted_doublets, num_predicted_doublets / num_total * 100.0],
    ]
    df = pd.DataFrame(data, columns=["value", "count", "% of total"])
    df["% of total"] = df["% of total"].round(2)
    st.dataframe(df)

    fig = get_scrublet_detected_doublet_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")


def get_scrublet_detected_doublet_distribution_plot(adata):
    fig, ax = plt.subplots()
    idx = adata.obs["predicted_doublet"]
    sns.histplot(
        adata.obs.loc[~idx, "doublet_score"], bins=50, stat="count", ax=ax, label="Singlets"
    )
    sns.histplot(
        adata.obs.loc[idx, "doublet_score"],
        bins=50,
        stat="count",
        color="red",
        ax=ax,
        label="Doublets",
    )
    ax.legend()
    ax.set_xlabel("Doublet Score")
    ax.set_ylabel("Barcode Count")
    ax.set_title("Detected Doublet Distribution")

    return fig


def display_highly_variable_genes_plot(adata):
    sc.pl.highly_variable_genes(adata)
    fig = plt.gcf()
    st.pyplot(fig)


def display_ranked_pca_importance(adata):
    # MRP
    data = [
        [n, variance_ratio] for n, variance_ratio in enumerate(adata.uns["pca"]["variance_ratio"])
    ]
    plot_df = pd.DataFrame(data, columns=["component", "variance_ratio"])

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax = sns.barplot(plot_df, x="component", y="variance_ratio", orient="x", ax=ax)
    ax.set_title("Ranked importance of Principal Components", fontsize=18)
    ax.set_ylabel("Variance Ratio", fontsize=20)
    ax.set_xlabel("Principal Component", fontsize=20)
    xticklabels = [
        "",
    ] * len(plot_df)
    for i in range(0, len(plot_df), 5):
        xticklabels[i] = str(i)
    ax.set_xticklabels(xticklabels)

    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def display_cumulative_pca_importance(adata):
    x = np.arange(0, len(adata.uns["pca"]["variance_ratio"]), 1)
    y = np.cumsum(adata.uns["pca"]["variance_ratio"])

    fig, ax = plt.subplots(figsize=(8, 4.75))

    ax.plot(x, y)
    ax.set_title("Cumulative importance of PCA components", fontsize=18)
    ax.set_xlabel("PCA component", fontsize=20)
    ax.set_ylabel("cumulative variance ratio", fontsize=20)

    # set a tick for every principal component, but only label every 5th
    ax.set_xticks(range(0, len(x)))
    xticklabels = [
        "",
    ] * len(x)
    for i in range(0, len(x), 5):
        xticklabels[i] = str(i)
    ax.set_xticklabels(xticklabels)

    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def get_top_principal_component_genes(adata, num_principal_components=3, keep_top_n=5):
    top_df_list = []
    bottom_df_list = []
    for pc_num in range(num_principal_components):
        # get indices that would sort principal component pc_num in descending order
        idx = np.argsort(adata.varm["PCs"][:, pc_num])
        idx = idx[-1::-1]

        # build a dataframe with the 'keep_top_n' largest principal component values
        top_genes = adata.var.iloc[idx[:keep_top_n]].index.tolist()
        top_values = adata.varm["PCs"][idx[:keep_top_n], pc_num].tolist()
        top_df = pd.DataFrame(data=zip(top_genes, top_values), columns=["gene", "PC value"])
        top_df_list.append(top_df)

        # build a dataframe with the 'keep_top_n' smallest principal component values
        bottom_genes = adata.var.iloc[idx[-keep_top_n - 1 :]].index.tolist()
        bottom_values = adata.varm["PCs"][idx[-keep_top_n - 1 :], pc_num].tolist()
        bottom_df = pd.DataFrame(
            data=zip(bottom_genes, bottom_values), columns=["gene", "PC value"]
        )
        bottom_df_list.append(bottom_df)

    # concatenate into two dataframes
    all_top_df = pd.concat(
        top_df_list, keys=[f"PC {i+1}" for i in range(num_principal_components)], axis=1
    )
    all_bottom_df = pd.concat(
        bottom_df_list, keys=[f"PC {i+1}" for i in range(num_principal_components)], axis=1
    )

    return all_top_df, all_bottom_df


def display_top_principal_component_genes(adata):
    all_top_df, all_bottom_df = get_top_principal_component_genes(
        adata, num_principal_components=3, keep_top_n=5
    )

    st.markdown("#### Genes with the Largest Principal Component Values")
    st.dataframe(all_top_df)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")

    st.markdown("#### Genes with the Smallest Principal Component Values")
    st.dataframe(all_bottom_df)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def display_projections_by_sample(adata):
    st.write("blah")

    sc.pl.pca(
        adata,
        color=["sample", "sample"],
        dimensions=[(0, 1), (2, 3)],
        size=2,
        title=["Sample Projection (PC1, PC2)", "Sample Projection (PC3, PC4)"],
    )
    fig = plt.gcf()
    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def display_projections_by_percent_mt(adata):
    sc.pl.pca(
        adata,
        color=["pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3)],
        size=2,
        title=["% MT Projection (PC1, PC2)", "% MT Projection (PC3, PC4)"],
    )
    fig = plt.gcf()
    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def display_umap_results(adata):
    sc.pl.umap(
        adata,
        color="sample",
        # Setting a smaller point size to get prevent overlap
        size=2,
    )
    fig = plt.gcf()
    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.markdown(
            "This plot provides a visualization of the entire dataset where each point corresponds to a barcode (i.e. a cell in most cases) and the colors correspond to the sample.  Barcodes that have similar expression profiles are near each other and barcodes that have different expression profiles appear far apart.  The exact shape can be tuned using the parameters to the left.  The algorithm is stochastic and so if you run with a different *random state*, you will get different results. This projection is based on the [McInnes et al., 2018](https://umap-learn.readthedocs.io/en/latest/) publication."
        )


def display_tsne_results(adata):
    sc.pl.tsne(
        adata,
        color="sample",
        # Setting a smaller point size to get prevent overlap
        size=2,
    )
    fig = plt.gcf()
    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


def display_clustering(adata, clustering_method, visualization_type):
    # set the label by which to color plots based on the clustering_method
    if clustering_method == "Leiden":
        label = "leiden"
    elif clustering_method == "Louvain":
        label = "louvain"
    else:
        raise ValueError(f"clustering_method '{clustering_method}' not recognized")

    # Plot the UMAP or t-SNE projection by clustering method
    if visualization_type == "UMAP":
        sc.pl.umap(adata, color=[label])
        ax = plt.gca()
        ax.set_title(f"{clustering_method} clustering", fontsize=20)
        fig = plt.gcf()
        st.pyplot(fig)

    elif visualization_type == "t-SNE":
        sc.pl.tsne(adata, color=[label])
        ax = plt.gca()
        ax = plt.gca()
        ax.set_title(f"{clustering_method} clustering", fontsize=20)
        fig = plt.gcf()
        st.pyplot(fig)

    else:
        raise ValueError(f"visualization_type '{visualization_type}' not recognized")

    with st.expander("More Info", expanded=False, icon="💭"):
        st.write(
            f"The clusters calculated using the {clustering_method} algorithm are shown on"
            f" the {visualization_type} projection calculated in the previous step.  Each"
            f" point corresponds to a barcode (i.e. a cell in most cases) and the colors"
            f" correspond to different clusters. Note that the cluster IDs (0, 1, ..) are"
            f" arbitrary. Barcodes within a cluster are generally"
            f" more similar to each other than barcodes from neighboring clusters. They"
            f" often signify groups with similar expression profiles and often distinguish"
            f" different cell types. Distances"
            f" between points should be interpreted with caution. See the caveats mentioned"
            f" in the previous step."
        )
