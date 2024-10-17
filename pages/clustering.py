from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import run_clustering
from utils.plotting import display_clustering


@dataclass
class ClusteringPageState:
    # parameters selected by user
    clustering_complete: bool = False
    clustering_method: Union[str, None] = None

    leiden_resolution: Union[float, None] = None
    leiden_n_iterations: Union[int, None] = None
    leiden_random_state: Union[int, None] = None

    louvain_resolution: Union[float, None] = None
    louvain_random_state: Union[int, None] = None

    # parameters not selected by user
    normalized_adata: Union[AnnData, None] = None
    visualization_type: Union[str, None] = None

    def reset(self):
        self.clustering_complete = False
        self.clustering_method = None
        self.leiden_resolution = None
        self.leiden_n_iterations = None
        self.leiden_random_state = None
        self.louvain_resolution = None
        self.louvain_random_state = None

    def save(self):
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


def display_clustering_method_options(clustering_method):
    page_state = st.session_state.clustering
    if clustering_method == "Leiden":

        page_state.user_sel_leiden_resolution = st.sidebar.number_input(
            "resolution",
            min_value=0.01,
            value=page_state.leiden_resolution if page_state.leiden_resolution else 1.0,
            help="A parameter that controls the coarseness of the clustering.  Larger values result"
            " in more clusters",
        )
        page_state.user_sel_leiden_n_iterations = st.sidebar.number_input(
            "Number of iterations",
            min_value=-1,
            value=page_state.leiden_n_iterations if page_state.leiden_n_iterations else 2,
            help="How many iterations of the algorithm to run.  -1 indicates it will run until it"
            " reaches an optimal clustering. (default: 2)",
        )
        page_state.user_sel_leiden_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=page_state.leiden_random_state if page_state.leiden_random_state else 0,
            help="Random state to use",
        )

    elif clustering_method == "Louvain":
        page_state.user_sel_louvain_resolution = st.sidebar.number_input(
            "resolution",
            min_value=0.01,
            value=page_state.louvain_resolution if page_state.louvain_resolution else 1.0,
            help="A parameter that controls the coarseness of the clustering.  Larger values result"
            " in more clusters",
        )
        page_state.user_sel_louvain_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=page_state.louvain_random_state if page_state.louvain_random_state else 0,
            help="Random state to use",
        )
    else:
        raise ValueError(f"cluster_method '{clustering_method}' not recognized")


def display_sidebar():
    # TODO: add KMeans using sklearn implementation
    page_state = st.session_state.clustering
    options_list = ["Leiden", "Louvain"]
    page_state.user_sel_clustering_method = st.sidebar.selectbox(
        "Clustering Method",
        options=options_list,
        index=(
            options_list.index(page_state.clustering_method) if page_state.clustering_method else 0
        ),
        help="Clustering Method to use.  The current recommendation of Scanpy and Seurat is to use"
        " the Leiden method."
        "\n\nLeiden: https://www.nature.com/articles/s41598-019-41695-z"
        "\n\nLouvain: https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008",
    )

    display_clustering_method_options(page_state.user_sel_clustering_method)

    run_clustering_clicked = st.sidebar.button("Run Clustering")

    if run_clustering_clicked:
        page_state.reset()
        page_state.save()
        run_clustering(page_state.normalized_adata, page_state)
        page_state.clustering_complete = True


def run():

    # add ClusteringPageState to track state if it hasn't already been initialized
    if "clustering" not in st.session_state:
        # if True:
        st.session_state.clustering = ClusteringPageState(
            normalized_adata=st.session_state.normalized_adata,
            visualization_type=st.session_state.visualization_type,
        )

    st.markdown("# Clustering")

    display_sidebar()

    if st.session_state.clustering.clustering_complete:
        display_clustering(
            st.session_state.clustering.normalized_adata,
            st.session_state.clustering.clustering_method,
            st.session_state.clustering.visualization_type,
        )
    else:
        st.write(
            "Clustering has not been run.  Select the desired parameters on the left and click Run"
            " Clustering."
        )


run()
