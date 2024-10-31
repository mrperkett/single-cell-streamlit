import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import run_dimensional_reduction_projection
from utils.plotting import display_tsne_results, display_umap_results
from utils.utils import set_downstream_pages_to_not_complete


@dataclass
class ProjectionPageState:
    # user selections when "Run Dimensional Reduction" button was last clicked
    visualization_type: Union[str, None] = None

    umap_min_dist: Union[float, None] = None
    umap_spread: Union[float, None] = None
    umap_random_state: Union[int, None] = None
    umap_alpha: Union[float, None] = None

    tsne_n_pcs: Union[int, None] = None
    tsne_perplexity: Union[float, None] = None
    tsne_early_exaggeration: Union[float, None] = None
    tsne_learning_rate: Union[int, None] = None
    tsne_random_state: Union[int, None] = None

    # current user selections
    user_sel_visualization_type: Union[str, None] = None

    user_sel_umap_min_dist: Union[float, None] = None
    user_sel_umap_spread: Union[float, None] = None
    user_sel_umap_random_state: Union[int, None] = None
    user_sel_umap_alpha: Union[float, None] = None

    user_sel_tsne_n_pcs: Union[int, None] = None
    user_sel_tsne_perplexity: Union[float, None] = None
    user_sel_tsne_early_exaggeration: Union[float, None] = None
    user_sel_tsne_learning_rate: Union[int, None] = None
    user_sel_tsne_random_state: Union[int, None] = None

    # parameters not selected by user
    normalized_adata: Union[AnnData, None] = None
    projection_complete: bool = False

    def reset(self):
        self.visualization_type = None
        self.umap_min_dist = None
        self.umap_spread = None
        self.umap_random_state = None
        self.umap_alpha = None
        self.tsne_n_pcs = None
        self.tsne_perplexity = None
        self.tsne_early_exaggeration = None
        self.tsne_learning_rate = None
        self.tsne_random_state = None
        self.projection_complete = False

    def update(self):
        self.visualization_type = self.user_sel_visualization_type

        if self.visualization_type == "UMAP":
            self.umap_min_dist = self.user_sel_umap_min_dist
            self.umap_spread = self.user_sel_umap_spread
            self.umap_random_state = self.user_sel_umap_random_state
            self.umap_alpha = self.user_sel_umap_alpha
        elif self.visualization_type == "t-SNE":
            self.tsne_n_pcs = self.user_sel_tsne_n_pcs
            self.tsne_perplexity = self.user_sel_tsne_perplexity
            self.tsne_early_exaggeration = self.user_sel_tsne_early_exaggeration
            self.tsne_learning_rate = self.user_sel_tsne_learning_rate
            self.tsne_random_state = self.user_sel_tsne_random_state
        else:
            raise ValueError(f"visualization type '{self.visualization_type}' not recognized")


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            try:
                self.state = ProjectionPageState(
                    normalized_adata=st.session_state.pca.normalized_adata
                )
            except Exception:
                self.state = ProjectionPageState()

    def display_umap_options(self):
        self.state.user_sel_umap_min_dist = st.sidebar.number_input(
            "Min distance",
            min_value=0.005,
            value=self.state.umap_min_dist if self.state.umap_min_dist else 0.05,
            help="The effective minimum distance between embedded points.  Smaller values result in"
            " a more tightly clustered embedding.  Larger values result in a more uniform"
            " distribution of points.",
        )
        self.state.user_sel_umap_spread = st.sidebar.number_input(
            "Spread",
            min_value=0.01,
            value=self.state.umap_spread if self.state.umap_spread else 1.0,
            help="The scale of embedded points.  Along with 'Min distance', this determines how"
            " tightly clustered the embedding is.",
        )
        self.state.user_sel_umap_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=self.state.umap_random_state if self.state.umap_random_state else 0,
            help="Random state to use",
        )
        self.state.user_sel_umap_alpha = st.sidebar.number_input(
            "Learning rate",
            min_value=0.01,
            value=self.state.umap_alpha if self.state.umap_alpha else 1.0,
            help="The initial learning rate for the embedding optimization (default: 1.0)",
        )

    def display_tsne_options(self):
        # Number of principal components was determined in the previous step.  Display it here,
        # but do not allow it to be modified
        self.state.user_sel_tsne_n_pcs = st.sidebar.number_input(
            "Number of Principal Components",
            min_value=1,
            # TODO: update value after PCA page has been refactored to contain a state class
            value=st.session_state.num_principal_components,
            help="Number of principal components to use use when running t-SNE (default: 50)",
            disabled=True,
        )

        self.state.user_sel_tsne_perplexity = st.sidebar.number_input(
            "Perplexity",
            min_value=5,
            max_value=50,
            value=self.state.tsne_perplexity if self.state.tsne_perplexity else 30,
            help="Related to the number of nearest neighbors. Larger datasets usually require"
            " larger perplexity (default: 30)",
        )

        self.state.user_sel_tsne_early_exaggeration = st.sidebar.number_input(
            "Early Exaggeration",
            min_value=0.01,
            value=(
                self.state.tsne_early_exaggeration if self.state.tsne_early_exaggeration else 12.0
            ),
            help="Controls how tight clusters are in the embedded space (default: 12)",
        )

        self.state.user_sel_tsne_learning_rate = st.sidebar.number_input(
            "Learning rate",
            min_value=100,
            value=self.state.tsne_learning_rate if self.state.tsne_learning_rate else 1000,
            help="The learning rate for the embedding optimization (default: 1000)",
        )

        self.state.user_sel_tsne_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=self.state.tsne_random_state if self.state.tsne_random_state else 0,
            help="Random state to use",
        )

    def display_sidebar(self):
        options_list = ["UMAP", "t-SNE"]
        self.state.user_sel_visualization_type = st.sidebar.radio(
            "Visualization Type",
            options=options_list,
            index=(
                options_list.index(self.state.visualization_type)
                if self.state.visualization_type
                else 0
            ),
            help="Select which projection to generate.",
        )

        if self.state.user_sel_visualization_type == "UMAP":
            self.display_umap_options()
        elif self.state.user_sel_visualization_type == "t-SNE":
            self.display_tsne_options()
        else:
            raise ValueError(
                f"visualization_type '{self.state.user_sel_visualization_type}' not " "recognized"
            )

        self.run_button_clicked = st.sidebar.button(
            "Run Dimensional Reduction",
            help="Run the dimensional reduction based on user parameters",
        )

    def run_dimensional_reduction_projection(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        # run UMAP / t-SNE projection
        with st.spinner("Running dimensional reduction projection"):
            run_dimensional_reduction_projection(self.state.normalized_adata, self.state)

        # save state to global st.session_state
        self.state.projection_complete = True
        self.save_to_session_state()

    def display_projection(self):
        if self.state.projection_complete:
            if self.state.visualization_type == "UMAP":
                display_umap_results(self.state.normalized_adata)
            elif self.state.visualization_type == "t-SNE":
                display_tsne_results(self.state.normalized_adata)
            else:
                raise ValueError(
                    f"visualization_type '{self.state.visualization_type}' not recognized"
                )
        else:
            st.markdown(
                "UMAP / t-SNE Projection has not been run.  Select the desired parameters on"
                " the left and click *Run Dimensional Reduction*"
            )

    def save_to_session_state(self):
        st.session_state.projection = self.state

    def run(self):
        st.markdown("# Dimensional Reduction: UMAP/t-SNE")
        page_step_number = st.session_state.page_completion_order.index("dim_reduction_umap_tsne")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()
        if self.run_button_clicked:
            self.run_dimensional_reduction_projection()
            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number
            set_downstream_pages_to_not_complete("dim_reduction_umap_tsne")
        self.display_projection()


if "projection" in st.session_state:
    page = Page(page_state=st.session_state.projection)
else:
    page = Page()
page.run()
