import matplotlib.pyplot as plt
import pandas as pd
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
    sns.histplot(adata.obs["total_counts"], ax=ax)
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
