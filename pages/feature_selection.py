import streamlit as st

from utils.analysis import mark_highly_variable_genes
from utils.plotting import display_highly_variable_genes_plot


def display_sidebar():
    st.session_state.feature_selection_algorithm = st.sidebar.selectbox(
        "Algorithm to detect highly variable genes",
        # options=["seurat", "cell_ranger", "seurat_v3"],
        options=["seurat", "cell_ranger"],  # drop seurat_v3 support for now, it needs tweaking
    )
    st.session_state.run_feature_selection_clicked = st.sidebar.button("Run Feature Selection")


def run():
    if "feature_selection_complete" not in st.session_state:
        st.session_state.feature_selection_complete = False

    st.markdown("# Feature Selection")

    display_sidebar()

    if st.session_state.run_feature_selection_clicked:
        mark_highly_variable_genes(
            st.session_state.normalized_adata, st.session_state.feature_selection_algorithm
        )
        st.session_state.feature_selection_complete = True

    if st.session_state.feature_selection_complete:
        display_highly_variable_genes_plot(st.session_state.normalized_adata)
    else:
        st.write("Feature selection has not been run.")


run()
