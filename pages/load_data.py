import time

import anndata as ad
import scanpy as sc
import streamlit as st
import yaml

from utils.plotting import (
    display_gene_distribution_plot,
    display_get_umi_distribution_plot,
    display_log_umi_count_log_gene_count_jointplot,
    display_mitochondrial_umi_distribution_plot,
    display_umi_count_gene_count_scatterplot,
)
from utils.preprocessing import run_quality_control
from utils.stats import get_qc_stats_dataframe


def load_config(config_file_path):
    with open(config_file_path, "r") as input_file:
        config = yaml.safe_load(input_file)
    return config


def load_dataset(dataset, status_placeholder):
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


def display_qc_info():

    # Display QC stats table
    df_general, df_umis_per_gene, df_umis_per_cell, df_genes_per_cell = get_qc_stats_dataframe(
        st.session_state.adata
    )
    st.dataframe(df_general, width=400)

    # (UMI count / cell) vs (gene count / cell) colored by (% MT)
    st.markdown("#")
    display_umi_count_gene_count_scatterplot(st.session_state.adata)

    # jointplot of (UMI count / cell) vs (gene count / cell)
    st.markdown("#")
    display_log_umi_count_log_gene_count_jointplot(st.session_state.adata)

    # Distribution of UMI count / cell
    st.markdown("#")
    display_get_umi_distribution_plot(st.session_state.adata)

    st.dataframe(df_umis_per_cell)

    # Distribution of gene count / cell
    st.markdown("#")
    display_gene_distribution_plot(st.session_state.adata)

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
    display_mitochondrial_umi_distribution_plot(st.session_state.adata)


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
