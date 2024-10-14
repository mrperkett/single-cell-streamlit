import numpy as np
import streamlit as st

from utils.plotting import (
    display_gene_distribution_plot,
    display_get_umi_distribution_plot,
    display_log_umi_count_log_gene_count_jointplot,
    display_mitochondrial_umi_distribution_plot,
    display_umi_count_gene_count_scatterplot,
)
from utils.stats import get_qc_stats_dataframe


def filter_adata(
    adata,
    min_allowed_genes_in_cell,
    min_allowed_cells_with_gene,
    max_allowed_percent_mt,
):
    # adata.obs["n_genes_by_counts"] gives total number of genes identified for each cell (i.e.
    # count all genes with >= 1 UMI)
    passes_min_allowed_genes_in_cell = adata.obs["n_genes_by_counts"] >= min_allowed_genes_in_cell

    passes_max_allowed_percent_mt = adata.obs["pct_counts_mt"] <= max_allowed_percent_mt

    # add barcode filters
    adata.obs["passes_filters"] = passes_min_allowed_genes_in_cell & passes_max_allowed_percent_mt

    # add feature filters
    # number of barcodes (i.e. cells) for each gene
    total_barcodes_by_gene = np.array((adata.X > 0).sum(axis=0)).flatten()
    adata.var["passes_filters"] = total_barcodes_by_gene >= min_allowed_cells_with_gene

    # create filtered_adata
    # filtered_adata = adata[adata.obs["passes_filters"]].copy()
    # filtered_adata = adata[:, adata.var["passes_filters"]].copy()
    filtered_adata = adata[adata.obs["passes_filters"], adata.var["passes_filters"]].copy()

    return filtered_adata


def define_sidebar_filtering_widgets():
    st.session_state.min_allowed_genes_in_cell = st.sidebar.number_input(
        "Min allowed genes in a cell",
        min_value=1,
        value="min",
        step=1,
        help="Cells that have less than this number of genes will be removed",
    )

    st.session_state.min_allowed_cells_with_gene = st.sidebar.number_input(
        "Min allowed cells with gene",
        min_value=1,
        value="min",
        step=1,
        help="Genes that are identified in less than this number of cells will be removed",
    )

    st.session_state.max_allowed_percent_mt = st.sidebar.slider(
        "Max allowed % Mitochondrial",
        min_value=0.0,
        max_value=100.0,
        value=100.0,
        step=1.0,
        help="Cells where the percent Mitochondrial UMIs is greater than this number will be removed",
    )

    st.session_state.apply_filters_button_clicked = st.sidebar.button(
        "Apply Filters",
        help="Apply quality control filters to data",
    )


def run():
    st.markdown("# Quality Control")

    # Define sidebar filtering selections
    define_sidebar_filtering_widgets()

    if "adata" not in st.session_state:
        st.markdown("**No data has been loaded.  Please run Load Data.**")
        return

    if "filtered_adata" not in st.session_state:
        st.session_state.filtered_adata = None

    # "Apply Filters" button clicked
    if st.session_state.apply_filters_button_clicked:
        # filter adata based on user selections
        st.session_state.filtered_adata = filter_adata(
            st.session_state.adata,
            st.session_state.min_allowed_genes_in_cell,
            st.session_state.min_allowed_cells_with_gene,
            st.session_state.max_allowed_percent_mt,
        )

    # load adata stats
    df_general, df_umis_per_gene, df_umis_per_cell, df_genes_per_cell = get_qc_stats_dataframe(
        st.session_state.adata
    )

    # load filtered_adata stats
    if st.session_state.filtered_adata:
        (
            df_general_filtered,
            df_umis_per_gene_filtered,
            df_umis_per_cell_filtered,
            df_genes_per_cell_filtered,
        ) = get_qc_stats_dataframe(st.session_state.filtered_adata)

    # define two columns for output

    # Display QC stats table
    column1, column2 = st.columns(2)
    with column1:
        st.markdown("## Original")
        st.dataframe(df_general, width=300)
    with column2:
        st.markdown("## Filtered")
        if st.session_state.filtered_adata:
            st.dataframe(df_general_filtered, width=300)
        else:
            st.markdown("Click *Apply Filters* to see filtered plots")

    # (UMI count / cell) vs (gene count / cell) colored by (% MT)
    st.markdown("")
    column1, column2 = st.columns(2)
    with column1:
        display_umi_count_gene_count_scatterplot(st.session_state.adata)
    with column2:
        if st.session_state.filtered_adata:
            display_umi_count_gene_count_scatterplot(st.session_state.filtered_adata)

    # jointplot of (UMI count / cell) vs (gene count / cell)
    st.markdown("")
    column1, column2 = st.columns(2)
    with column1:
        display_log_umi_count_log_gene_count_jointplot(st.session_state.adata)
    with column2:
        if st.session_state.filtered_adata:
            display_log_umi_count_log_gene_count_jointplot(st.session_state.filtered_adata)

    # Distribution of UMI count / cell
    st.markdown("#")
    column1, column2 = st.columns(2)
    with column1:
        display_get_umi_distribution_plot(st.session_state.adata)
        st.dataframe(df_umis_per_cell)
    with column2:
        if st.session_state.filtered_adata:
            display_get_umi_distribution_plot(st.session_state.filtered_adata)
            st.dataframe(df_umis_per_cell_filtered)

    # Distribution of gene count / cell
    st.markdown("#")
    column1, column2 = st.columns(2)
    with column1:
        display_gene_distribution_plot(st.session_state.adata)
        st.dataframe(df_genes_per_cell)
    with column2:
        if st.session_state.filtered_adata:
            display_gene_distribution_plot(st.session_state.filtered_adata)
            st.dataframe(df_genes_per_cell_filtered)

    # Distribution of % MT UMIs
    st.markdown("#")
    column1, column2 = st.columns(2)
    with column1:
        display_mitochondrial_umi_distribution_plot(st.session_state.adata)
    with column2:
        if st.session_state.filtered_adata:
            display_mitochondrial_umi_distribution_plot(st.session_state.filtered_adata)


run()
