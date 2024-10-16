import streamlit as st

from utils.analysis import run_pca
from utils.plotting import (
    display_cumulative_pca_importance,
    display_projections_by_percent_mt,
    display_projections_by_sample,
    display_ranked_pca_importance,
    display_top_principal_component_genes,
)


def display_sidebar():
    st.session_state.num_principal_components = st.sidebar.number_input(
        "Num Principal Components",
        min_value=1,
        value=50,
        step=1,
        help="The number of principal components to keep",
    )
    st.session_state.pca_use_highly_variable = st.sidebar.checkbox(
        "Use Highly Variable Genes",
        value=True,
        help="Only use the highly variable genes identified in the 'Feature Selection' step",
    )
    st.session_state.pca_random_state = st.sidebar.number_input(
        "Random State", min_value=0, value=0, step=1, help="The random seed used"
    )

    run_pca_button_clicked = st.sidebar.button(
        "Run PCA",
        help="Run principal component analysis based on user selections and generate plots.",
    )

    if run_pca_button_clicked:
        run_pca(
            st.session_state.normalized_adata,
            n_comps=st.session_state.num_principal_components,
            use_highly_variable=st.session_state.pca_use_highly_variable,
            random_state=st.session_state.pca_random_state,
        )
        st.session_state.run_pca_complete = True


def run():
    if "run_pca_complete" not in st.session_state:
        st.session_state.run_pca_complete = False

    st.markdown("# Dimensional Reduction")

    display_sidebar()

    if st.session_state.run_pca_complete:
        display_ranked_pca_importance(st.session_state.normalized_adata)

        st.markdown("##")
        display_cumulative_pca_importance(st.session_state.normalized_adata)

        st.markdown("##")
        display_top_principal_component_genes(st.session_state.normalized_adata)

        st.markdown("#### Principal Component Projections")
        display_projections_by_sample(st.session_state.normalized_adata)
        display_projections_by_percent_mt(st.session_state.normalized_adata)
    else:
        st.write("PCA has not been run.  Click 'Run PCA' to continue.")


run()
