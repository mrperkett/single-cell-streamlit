import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import mark_highly_variable_genes
from utils.plotting import display_highly_variable_genes_plot


@dataclass
class FeatureSelectionPageState:
    # user selections when "Run Dimensional Reduction" button was last clicked
    feature_selection_algorithm: Union[str, None] = None

    # current user selections
    user_sel_feature_selection_algorithm: Union[str, None] = None

    # parameters not selected by user
    normalized_adata: Union[AnnData, None] = None
    feature_selection_complete: bool = False

    def reset(self):
        self.feature_selection_algorithm = None

    def update(self):
        self.feature_selection_algorithm = self.user_sel_feature_selection_algorithm


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            # TODO: update normalized_adata location when previous steps store this as
            # st.session_state.normalization.normalized_adata
            self.state = FeatureSelectionPageState(
                normalized_adata=st.session_state.normalized_adata
            )

    def display_sidebar(self):
        # TODO: drop seurat_v3 support for now since it needs tweaking
        options_list = ["seurat", "cell_ranger"]
        self.state.user_sel_feature_selection_algorithm = st.sidebar.selectbox(
            "Algorithm to detect highly variable genes",
            options=options_list,
            index=(
                options_list.index(self.state.feature_selection_algorithm)
                if self.state.feature_selection_algorithm
                else 0
            ),
        )
        self.run_feature_selection_clicked = st.sidebar.button("Run Feature Selection")

    def run_feature_selection(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        mark_highly_variable_genes(
            st.session_state.normalized_adata, st.session_state.feature_selection_algorithm
        )

        # save state to global st.session_state
        st.session_state.feature_selection_complete = True
        self.save_to_session_state()

    def save_to_session_state(self):
        st.session_state.feature_selection = self.state

    def run(self):
        st.markdown("# Feature Selection")
        page_step_number = st.session_state.page_completion_order.index("dim_reduction_umap_tsne")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()

        if self.run_feature_selection_clicked:
            self.run_feature_selection()
            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number
            if "pca" in st.session_state:
                st.session_state.pca.run_pca_complete = False
            if "projection" in st.session_state:
                st.session_state.projection.projection_complete = False
            if "clustering" in st.session_state:
                st.session_state.clustering.clustering_complete = False

        if st.session_state.feature_selection_complete:
            display_highly_variable_genes_plot(st.session_state.normalized_adata)
        else:
            st.write("Feature selection has not been run.")


if "feature_selection" in st.session_state:
    page = Page(page_state=st.session_state.feature_selection)
else:
    page = Page()
page.run()
