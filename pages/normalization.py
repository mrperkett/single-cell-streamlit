import streamlit as st

from utils.analysis import run_normalization
from utils.plotting import display_get_umi_distribution_plot


def display_sidebar():
    st.session_state.exclude_highly_expressed = st.sidebar.checkbox(
        "Exclude Highly Expressed Genes",
        value=True,
        help="From Scanpy documentation: "
        "Exclude (very) highly expressed genes for the computation of the normalization factor"
        " (size factor) for each cell. A gene is considered highly expressed, if it has more "
        "than max_fraction of the total counts in at least one cell. The not-excluded genes"
        "will sum up to target_sum. Providing this argument when adata.X is a Array will incur"
        " blocking .compute() calls on the array. [From scanpy documentation]",
    )
    st.session_state.max_fraction = st.sidebar.number_input(
        "Max Fraction",
        min_value=0.0,
        max_value=1.0,
        value=0.05,
        step=0.01,
        help="From Scanpy documentation: "
        "If exclude_highly_expressed=True, consider cells as highly expressed that have more "
        "counts than max_fraction of the original total counts in at least one cell.",
        disabled=not st.session_state.exclude_highly_expressed,
    )
    st.session_state.target_sum_category = st.sidebar.selectbox(
        "Target Sum",
        options=["Median", "Custom"],
        help="From Scanpy documentation: "
        "If None, after normalization, each observation (cell) has a total count equal to the"
        " median of total counts for observations (cells) before normalization.",
    )

    st.session_state.target_sum = st.sidebar.number_input(
        "Target Sum Custom Value",
        value=10000,
        disabled=st.session_state.target_sum_category != "Custom",
    )

    st.session_state.run_normalization_button_clicked = st.sidebar.button(
        "Run Normalization",
    )


def run():

    if "run_normalization_complete" not in st.session_state:
        st.session_state.run_normalization_complete = False

    st.markdown("# Normalization")

    # Sidebar options
    display_sidebar()

    # Run normalization
    if st.session_state.run_normalization_button_clicked:
        if st.session_state.target_sum_category == "Median":
            target_sum = None
        else:
            target_sum = st.session_state.target_sum
        st.session_state.normalized_adata = run_normalization(
            st.session_state.doublets_removed_adata,
            st.session_state.exclude_highly_expressed,
            st.session_state.max_fraction,
            target_sum,
        )

        st.session_state.run_normalization_complete = True

    if st.session_state.run_normalization_complete:
        st.write("Data has been normalized!")
        col1, col2 = st.columns(2)

        with col1:
            st.write("## Unnormalized")
            display_get_umi_distribution_plot(st.session_state.doublets_removed_adata)
        with col2:
            st.write("## Normalized")
            display_get_umi_distribution_plot(st.session_state.normalized_adata)


run()
