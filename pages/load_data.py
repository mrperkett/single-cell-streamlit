import time

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import streamlit as st
import yaml


def load_config(config_file_path):
    with open(config_file_path, "r") as input_file:
        config = yaml.safe_load(input_file)
    return config


def load_dataset(dataset, status_placeholder):
    # progress_bar = st.progress(0, text="Loading data, please wait.")
    # progress_bar.progress(float(n / num_samples), f"Loading sample {sample_id}")

    with status_placeholder.container():
        with st.status("Loading data", expanded=True) as status:
            adatas = {}
            # num_samples = len(dataset["samples"])
            for n, (sample_id, sample_path) in enumerate(dataset["samples"].items()):
                st.write(f"Loading sample {sample_id}")
                sample_adata = sc.read_10x_h5(sample_path)
                sample_adata.var_names_make_unique()
                adatas[sample_id] = sample_adata
                time.sleep(1)

            st.write("Concatenating samples")
            adata = ad.concat(adatas, label="sample")
            adata.obs_names_make_unique()
            time.sleep(1)

            st.session_state.adata = adata
            st.session_state.dataset_loaded = True
            st.session_state.quality_control_complete = False

        status.update(label="Load complete!", state="complete", expanded=False)


def run_quality_control(adata):
    # annotate genes: mitochondrial ("mt"), ribosomal ("rb"), hemoglobin ("hb")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    # run QC
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    st.session_state.quality_control_complete = True


def get_qc_stats_dataframe(adata):
    # total UMIs for each gene
    total_umis_by_gene = adata.var["total_counts"].astype(int)
    idx = total_umis_by_gene > 0
    min_umis_per_gene = np.min(total_umis_by_gene[idx])
    median_umis_per_gene = np.median(total_umis_by_gene[idx])
    mean_umis_per_gene = np.mean(total_umis_by_gene[idx])
    max_umis_per_gene = np.max(total_umis_by_gene[idx])
    stdev_umis_per_gene = np.std(total_umis_by_gene[idx])

    # total UMIs for each cell
    total_umis_by_cell = adata.obs["total_counts"].astype(int)
    idx = total_umis_by_cell > 0
    min_umis_per_cell = np.min(total_umis_by_cell[idx])
    median_umis_per_cell = np.median(total_umis_by_cell[idx])
    mean_umis_per_cell = np.mean(total_umis_by_cell[idx])
    max_umis_per_cell = np.max(total_umis_by_cell[idx])
    stdev_umis_per_cell = np.std(total_umis_by_cell[idx])

    # total genes for each cell
    total_genes_by_cell = adata.obs["n_genes_by_counts"]
    idx = total_genes_by_cell > 0
    min_genes_per_cell = np.min(total_genes_by_cell[idx])
    median_genes_per_cell = np.median(total_genes_by_cell[idx])
    mean_genes_per_cell = np.mean(total_genes_by_cell[idx])
    max_genes_per_cell = np.max(total_genes_by_cell[idx])
    stdev_genes_per_cell = np.std(total_genes_by_cell[idx])

    # totals over
    total_umi_count = adata.X.astype(int).sum()
    total_genes_identified = np.count_nonzero(total_umis_by_gene)

    data = [
        ["Total UMI Count", total_umi_count],
        ["Total Genes Identified", total_genes_identified],
        ["Median UMIs per Gene", median_umis_per_gene],
        ["Mean UMIs per Gene", mean_umis_per_gene],
        ["Median UMIs per Cell", median_umis_per_cell],
        ["Mean UMIs per Cell", mean_umis_per_cell],
    ]
    df_general = pd.DataFrame(data, columns=["Label", "Value"])

    data = [
        ["Min UMIs per Gene", min_umis_per_gene],
        ["Median UMIs per Gene", median_umis_per_gene],
        ["Mean UMIs per Gene", mean_umis_per_gene],
        ["Max UMIs per Gene", max_umis_per_gene],
        ["Stev UMIs per Gene", stdev_umis_per_gene],
    ]
    df_umis_per_gene = pd.DataFrame(data, columns=["Label", "Value"])

    data = [
        ["Min UMIs per Cell", min_umis_per_cell],
        ["Median UMIs per Cell", median_umis_per_cell],
        ["Mean UMIs per Cell", mean_umis_per_cell],
        ["Max UMIs per Cell", max_umis_per_cell],
        ["Stdev UMIs per Cell", stdev_umis_per_cell],
    ]
    df_umis_per_cell = pd.DataFrame(data, columns=["Label", "Value"])

    data = [
        ["Min genes per Cell", min_genes_per_cell],
        ["Median genes per Cell", median_genes_per_cell],
        ["Mean genes per Cell", mean_genes_per_cell],
        ["Max genes per Cell", max_genes_per_cell],
        ["Stdev genes per Cell", stdev_genes_per_cell],
    ]
    df_genes_per_cell = pd.DataFrame(data, columns=["Label", "Value"])

    return df_general, df_umis_per_gene, df_umis_per_cell, df_genes_per_cell


def display_qc_info():

    adata = st.session_state.adata

    # Display QC stats table
    df_general, df_umis_per_gene, df_umis_per_cell, df_genes_per_cell = get_qc_stats_dataframe(
        adata
    )
    st.dataframe(df_general, width=400)

    # (UMI count / cell) vs (gene count / cell) colored by (% MT)
    st.markdown("#")
    fig = get_umi_count_gene_count_scatterplot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")

    # jointplot of (UMI count / cell) vs (gene count / cell)
    st.markdown("#")
    fig = get_log_umi_count_log_gene_count_jointplot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")

    # Distribution of UMI count / cell
    st.markdown("#")
    fig = get_umi_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")

    st.dataframe(df_umis_per_cell)

    # Distribution of gene count / cell
    st.markdown("#")
    fig = get_gene_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander(" More Info", expanded=False, icon="💭"):
        st.write("description")

    st.dataframe(df_genes_per_cell)

    # NOTE: this runs really slow - exclude for now
    # # Distribution of UMIs / gene
    # st.markdown("#")
    # fig = get_umis_per_gene_distribution_plot(adata)
    # st.pyplot(fig)
    # with st.expander("More Info", expanded=False, icon="💭"):
    #     st.write("description")

    # st.dataframe(df_umis_per_gene)

    # Distribution of % MT UMIs
    st.markdown("#")
    fig = get_mitochondrial_umi_distribution_plot(adata)
    st.pyplot(fig)
    with st.expander("More Info", expanded=False, icon="💭"):
        st.write("description")


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


# def get_umis_per_gene_distribution_plot(adata):
#     fig, ax = plt.subplots()

#     sns.histplot(adata.var["total_counts"].to_numpy(), ax=ax)
#     ax.set_xlabel("UMI Count", fontsize=16)
#     ax.set_ylabel("Gene Count", fontsize=16)
#     ax.set_title("Distribution for UMIs / Gene", fontsize=18)

#     fig = plt.gcf()
#     return fig


def get_mitochondrial_umi_distribution_plot(adata):
    fig, ax = plt.subplots()
    sns.histplot(adata.obs["pct_counts_mt"], ax=ax)
    ax.set_xlabel("% UMIs that are Mitochondrial", fontsize=16)
    ax.set_ylabel("Barcode Count", fontsize=16)
    ax.set_title("Distribution of % Mitochondrial", fontsize=18)

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


def run():
    if "dataset_loaded" not in st.session_state:
        st.session_state.dataset_loaded = False
    if "quality_control_complete" not in st.session_state:
        st.session_state.quality_control_complete = False

    st.markdown("# Load Data")

    config = load_config("config.yaml")

    st.markdown("## Select Dataset")
    dataset_names = config["datasets"].keys()
    dataset_name = st.selectbox("Dataset", options=dataset_names)
    dataset = config["datasets"][dataset_name]

    description = f"*{dataset['description']}*"
    st.markdown(f"{description}")

    # Run data loading on button click
    load_button_clicked = st.button("Load")
    status_placeholder = st.empty()
    if load_button_clicked:
        st.session_state.dataset_name = dataset_name
        st.session_state.dataset = dataset
        load_dataset(dataset, status_placeholder)

    if st.session_state.dataset_loaded:
        st.write(f"Dataset '{st.session_state.dataset_name}' is loaded.")

    st.markdown("## Quality Control")

    # Run Quality Control if the dataset has been loaded and it hasn't already been run
    if st.session_state.dataset_loaded and not st.session_state.quality_control_complete:
        with st.spinner("Running QC..."):
            time.sleep(1)
            run_quality_control(st.session_state.adata)

    # Display Quality Control information after it has been run
    if st.session_state.quality_control_complete:
        display_qc_info()


run()
