import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import run_clustering
from utils.plotting import display_clustering
from utils.utils import set_downstream_pages_to_not_complete


@dataclass
class ClusteringPageState:
    # user selections when clustering was last run
    clustering_method: Union[str, None] = None

    leiden_resolution: Union[float, None] = None
    leiden_n_iterations: Union[int, None] = None
    leiden_random_state: Union[int, None] = None

    louvain_resolution: Union[float, None] = None
    louvain_random_state: Union[int, None] = None

    # current user selections
    user_sel_clustering_method: Union[str, None] = None

    user_sel_leiden_resolution: Union[float, None] = None
    user_sel_leiden_n_iterations: Union[int, None] = None
    user_sel_leiden_random_state: Union[int, None] = None

    user_sel_louvain_resolution: Union[float, None] = None
    user_sel_louvain_random_state: Union[int, None] = None

    # parameters not selected by user
    clustering_complete: bool = False
    normalized_adata: Union[AnnData, None] = None
    visualization_type: Union[str, None] = None

    def reset(self):
        self.clustering_method = None
        self.leiden_resolution = None
        self.leiden_n_iterations = None
        self.leiden_random_state = None
        self.louvain_resolution = None
        self.louvain_random_state = None
        self.clustering_complete = False

    def update(self):
        self.clustering_method = self.user_sel_clustering_method

        if self.clustering_method == "Leiden":
            self.leiden_resolution = self.user_sel_leiden_resolution
            self.leiden_n_iterations = self.user_sel_leiden_n_iterations
            self.leiden_random_state = self.user_sel_leiden_random_state
        elif self.clustering_method == "Louvain":
            self.louvain_resolution = self.user_sel_louvain_resolution
            self.louvain_random_state = self.user_sel_louvain_random_state
        else:
            raise ValueError(f"clustering method '{self.clustering_method}' not recognized")


class Page:
    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            try:
                self.state = ClusteringPageState(
                    normalized_adata=st.session_state.projection.normalized_adata,
                    visualization_type=st.session_state.projection.visualization_type,
                )
            except Exception:
                self.state = ClusteringPageState()
        self.run_clustering_clicked = False

    def display_leiden_options(self):
        self.state.user_sel_leiden_resolution = st.sidebar.number_input(
            "resolution",
            min_value=0.01,
            value=self.state.leiden_resolution if self.state.leiden_resolution else 1.0,
            help="A parameter that controls the coarseness of the clustering.  Larger values result"
            " in more clusters",
        )
        self.state.user_sel_leiden_n_iterations = st.sidebar.number_input(
            "Number of iterations",
            min_value=-1,
            value=self.state.leiden_n_iterations if self.state.leiden_n_iterations else 2,
            help="How many iterations of the algorithm to run.  -1 indicates it will run until it"
            " reaches an optimal clustering. (default: 2)",
        )
        self.state.user_sel_leiden_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=self.state.leiden_random_state if self.state.leiden_random_state else 0,
            help="Random state to use",
        )

    def display_louvain_options(self):
        self.state.user_sel_louvain_resolution = st.sidebar.number_input(
            "resolution",
            min_value=0.01,
            value=self.state.louvain_resolution if self.state.louvain_resolution else 1.0,
            help="A parameter that controls the coarseness of the clustering.  Larger values result"
            " in more clusters",
        )
        self.state.user_sel_louvain_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=self.state.louvain_random_state if self.state.louvain_random_state else 0,
            help="Random state to use",
        )

    def display_clustering_method_options(self, clustering_method):
        if clustering_method == "Leiden":
            self.display_leiden_options()
        elif clustering_method == "Louvain":
            self.display_louvain_options()
        else:
            raise ValueError(f"cluster_method '{clustering_method}' not recognized")

    def display_sidebar(self):
        # TODO: add KMeans using sklearn implementation
        options_list = ["Leiden", "Louvain"]
        self.state.user_sel_clustering_method = st.sidebar.selectbox(
            "Clustering Method",
            options=options_list,
            index=(
                options_list.index(self.state.clustering_method)
                if self.state.clustering_method
                else 0
            ),
            help="Clustering Method to use.  The current recommendation of Scanpy and Seurat is to"
            " use the Leiden method."
            "\n\nLeiden: https://www.nature.com/articles/s41598-019-41695-z"
            "\n\nLouvain: https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008",
        )

        self.display_clustering_method_options(self.state.user_sel_clustering_method)

        self.run_clustering_clicked = st.sidebar.button("Run Clustering")

    def save_to_session_state(self):
        st.session_state.clustering = self.state

    def run_clustering(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        # run clustering
        run_clustering(self.state.normalized_adata, self.state)

        # save state to global st.session_state
        self.state.clustering_complete = True
        self.save_to_session_state()

    def display_clustering(self):
        if self.state.clustering_complete:
            display_clustering(
                self.state.normalized_adata,
                self.state.clustering_method,
                self.state.visualization_type,
            )
        else:
            st.write(
                "Clustering has not been run.  Select the desired parameters on the left and click"
                " *Run Clustering*."
            )

    def run(self):
        st.markdown("# Clustering")
        page_step_number = st.session_state.page_completion_order.index("clustering")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()
        if self.run_clustering_clicked:
            self.run_clustering()
            st.session_state.furthest_step_number_completed = page_step_number
            set_downstream_pages_to_not_complete("clustering")
        self.display_clustering()


if "clustering" in st.session_state:
    page = Page(page_state=st.session_state.clustering)
else:
    page = Page()
page.run()
