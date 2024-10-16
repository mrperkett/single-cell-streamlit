import streamlit as st

from utils.analysis import run_clustering
from utils.plotting import display_clustering


def display_clustering_method_options(clustering_method):

    if clustering_method == "Leiden":

        st.session_state.current_leiden_resolution = st.sidebar.number_input(
            "resolution",
            min_value=0.01,
            value=1.0,
            help="A parameter that controls the coarseness of the clustering.  Larger values result"
            " in more clusters",
        )
        st.session_state.current_leiden_n_iterations = st.sidebar.number_input(
            "Number of iterations",
            min_value=-1,
            value=2,
            help="How many iterations of the algorithm to run.  -1 indicates it will run until it"
            " reaches an optimal clustering. (default: 2)",
        )
        st.session_state.current_leiden_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=0,
            help="Random state to use",
        )

    elif clustering_method == "Louvain":
        st.session_state.current_louvain_resolution = st.sidebar.number_input(
            "resolution",
            min_value=0.01,
            value=1.0,
            help="A parameter that controls the coarseness of the clustering.  Larger values result"
            " in more clusters",
        )
        st.session_state.current_louvain_random_state = st.sidebar.number_input(
            "Random state",
            min_value=0,
            value=0,
            help="Random state to use",
        )
    else:
        raise ValueError(f"cluster_method '{clustering_method}' not recognized")


def display_sidebar():
    # TODO: add KMeans using sklearn implementation
    current_clustering_method = st.sidebar.selectbox(
        "Clustering Method",
        options=["Leiden", "Louvain"],
        help="Clustering Method to use.  The current recommendation of Scanpy and Seurat is to use"
        " the Leiden method."
        "\n\nLeiden: https://www.nature.com/articles/s41598-019-41695-z"
        "\n\nLouvain: https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008",
    )

    display_clustering_method_options(current_clustering_method)

    run_clustering_clicked = st.sidebar.button("Run Clustering")

    if run_clustering_clicked:
        st.session_state.clustering_method = current_clustering_method
        run_clustering(st.session_state.normalized_adata, st.session_state)
        st.session_state.clustering_complete = True


def run():

    if "clustering_complete" not in st.session_state:
        st.session_state.clustering_complete = False

    st.markdown("# Clustering")

    display_sidebar()

    if st.session_state.clustering_complete:
        display_clustering(
            st.session_state.normalized_adata,
            st.session_state.clustering_method,
            st.session_state.visualization_type,
        )
    else:
        st.write(
            "Clustering has not been run.  Select the desired parameters on the left and click Run"
            " Clustering."
        )


run()
