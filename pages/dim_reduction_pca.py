import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import run_pca
from utils.plotting import (
    display_cumulative_pca_importance,
    display_projections_by_percent_mt,
    display_projections_by_sample,
    display_ranked_pca_importance,
    display_top_principal_component_genes,
)


@dataclass
class PCAPageState:
    # user selections when "Run Dimensional Reduction" button was last clicked
    num_principal_components: Union[int, None] = None
    pca_use_highly_variable: Union[bool, None] = None
    pca_random_state: Union[int, None] = None

    # current user selections
    user_sel_num_principal_components: Union[int, None] = None
    user_sel_pca_use_highly_variable: Union[bool, None] = None
    user_sel_pca_random_state: Union[int, None] = None

    # parameters not selected by user
    normalized_adata: Union[AnnData, None] = None
    run_pca_complete: bool = False

    def reset(self):
        self.num_principal_components = None
        self.pca_use_highly_variable = None
        self.pca_random_state = None
        self.run_pca_complete = False

    def update(self):
        self.num_principal_components = self.user_sel_num_principal_components
        self.pca_use_highly_variable = self.user_sel_pca_use_highly_variable
        self.pca_random_state = self.user_sel_pca_random_state


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            # TODO: update normalized_adata location when previous steps store this as
            # st.session_state.normalization.normalized_adata
            self.state = PCAPageState(normalized_adata=st.session_state.normalized_adata)
        self.run_pca_button_clicked = False

    def display_sidebar(self):
        self.state.user_sel_num_principal_components = st.sidebar.number_input(
            "Num Principal Components",
            min_value=1,
            value=(
                self.state.num_principal_components if self.state.num_principal_components else 50
            ),
            step=1,
            help="The number of principal components to keep",
        )
        self.state.user_sel_pca_use_highly_variable = st.sidebar.checkbox(
            "Use Highly Variable Genes",
            value=(
                self.state.pca_use_highly_variable if self.state.pca_use_highly_variable else True
            ),
            help="Only use the highly variable genes identified in the 'Feature Selection' step",
        )
        self.state.user_sel_pca_random_state = st.sidebar.number_input(
            "Random State",
            min_value=0,
            value=self.state.pca_random_state if self.state.pca_random_state else 0,
            step=1,
            help="The random seed used",
        )

        self.run_pca_button_clicked = st.sidebar.button(
            "Run PCA",
            help="Run principal component analysis based on user selections and generate plots.",
        )

    def display_pca_results(self):
        if self.state.run_pca_complete:
            display_ranked_pca_importance(self.state.normalized_adata)

            st.markdown("##")
            display_cumulative_pca_importance(self.state.normalized_adata)

            st.markdown("##")
            display_top_principal_component_genes(self.state.normalized_adata)

            st.markdown("#### Principal Component Projections")
            display_projections_by_sample(self.state.normalized_adata)
            display_projections_by_percent_mt(self.state.normalized_adata)
        else:
            st.markdown(
                "PCA has not been run.  Select the desired parameters on the left and"
                " click *Run PCA* to continue."
            )

    def run_pca(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        run_pca(
            self.state.normalized_adata,
            n_comps=self.state.num_principal_components,
            use_highly_variable=self.state.pca_use_highly_variable,
            random_state=self.state.pca_random_state,
        )

        # save state to global st.session_state
        self.state.run_pca_complete = True
        self.save_to_session_state()

    def save_to_session_state(self):
        st.session_state.pca = self.state

    def run(self):
        st.markdown("# Dimensional Reduction: PCA")
        page_step_number = st.session_state.page_completion_order.index("dim_reduction_pca")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()

        if self.run_pca_button_clicked:
            self.run_pca()

            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number
            if "projection" in st.session_state:
                st.session_state.projection.projection_complete = False
            if "clustering" in st.session_state:
                st.session_state.clustering.clustering_complete = False

        self.display_pca_results()


if "pca" in st.session_state:
    page = Page(page_state=st.session_state.pca)
else:
    page = Page()
page.run()
