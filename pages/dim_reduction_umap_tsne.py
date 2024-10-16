import streamlit as st

from utils.analysis import run_tsne, run_umap
from utils.plotting import display_tsne_results, display_umap_results


def display_umap_options():
    st.session_state.current_umap_min_dist = st.sidebar.number_input(
        "Min distance",
        min_value=0.005,
        value=0.05,
        help="The effective minimum distance between embedded points.  Smaller values result in a"
        " more tightly clustered embedding.  Larger values result in a more uniform distribution"
        " of points.",
    )
    st.session_state.current_umap_spread = st.sidebar.number_input(
        "Spread",
        min_value=0.01,
        value=1.0,
        help="The scale of embedded points.  Along with 'Min distance', this determines how tightly"
        " clustered the embedding is.",
    )
    st.session_state.current_umap_random_state = st.sidebar.number_input(
        "Random state",
        min_value=0,
        value=0,
        help="Random state to use",
    )
    st.session_state.current_umap_alpha = st.sidebar.number_input(
        "Learning rate",
        min_value=0.01,
        value=1.0,
        help="The initial learning rate for the embedding optimization (default: 1.0)",
    )


def display_tsne_options():
    # Number of principal components was determined in the previous step.  Display it here,
    # but do not allow it to be modified
    st.session_state.current_tsne_n_pcs = st.sidebar.number_input(
        "Number of Principal Components",
        min_value=1,
        value=st.session_state.num_principal_components,
        help="Number of principal components to use use when running t-SNE (default: 50)",
        disabled=True,
    )

    st.session_state.current_tsne_perplexity = st.sidebar.number_input(
        "Perplexity",
        min_value=5,
        max_value=50,
        value=30,
        help="Related to the number of nearest neighbors. Larger datasets usually require larger"
        " perplexity (default: 30)",
    )

    st.session_state.current_tsne_early_exaggeration = st.sidebar.number_input(
        "Early Exaggeration",
        min_value=0.01,
        value=12.0,
        help="Controls how tight clusters are in the embedded space (default: 12)",
    )

    st.session_state.current_tsne_learning_rate = st.sidebar.number_input(
        "Learning rate",
        min_value=100,
        value=1000,
        help="The learning rate for the embedding optimization (default: 1000)",
    )

    st.session_state.current_tsne_random_state = st.sidebar.number_input(
        "Random state",
        min_value=0,
        value=0,
        help="Random state to use",
    )


def display_sidebar():
    st.session_state.current_visualization_type = st.sidebar.radio(
        "Visualization Type",
        options=["UMAP", "t-SNE"],
        help="",
    )

    if st.session_state.current_visualization_type == "UMAP":
        display_umap_options()
    elif st.session_state.current_visualization_type == "t-SNE":
        display_tsne_options()
    else:
        raise ValueError(
            f"visualization_type '{st.session_state.current_visualization_type}' not recognized"
        )

    run_dim_red_visualization = st.sidebar.button(
        "Run Dimensional Reduction",
        help="Run the dimensional reduction based on user parameters",
    )

    if run_dim_red_visualization:
        st.session_state.visualization_type = st.session_state.current_visualization_type
        if st.session_state.visualization_type == "UMAP":
            # store current user selections that will be used for the run
            st.session_state.umap_min_dist = st.session_state.current_umap_min_dist
            st.session_state.umap_spread = st.session_state.current_umap_spread
            st.session_state.umap_random_state = st.session_state.current_umap_random_state
            st.session_state.umap_alpha = st.session_state.current_umap_alpha

            # run UMAP
            with st.spinner("Running dimensional reduction"):
                run_umap(
                    adata=st.session_state.normalized_adata,
                    min_dist=st.session_state.umap_min_dist,
                    spread=st.session_state.umap_spread,
                    random_state=st.session_state.umap_random_state,
                    alpha=st.session_state.umap_alpha,
                )
        elif st.session_state.visualization_type == "t-SNE":

            # store current user selections that will be used for the run
            st.session_state.tsne_n_pcs = st.session_state.current_tsne_n_pcs
            st.session_state.tsne_perplexity = st.session_state.current_tsne_perplexity
            st.session_state.tsne_early_exaggeration = (
                st.session_state.current_tsne_early_exaggeration
            )
            st.session_state.tsne_learning_rate = st.session_state.current_tsne_learning_rate
            st.session_state.tsne_random_state = st.session_state.current_tsne_random_state

            # run t-SNE
            with st.spinner(
                "Running dimensional reduction",
            ):
                run_tsne(
                    adata=st.session_state.normalized_adata,
                    n_pcs=st.session_state.tsne_n_pcs,
                    perplexity=st.session_state.tsne_perplexity,
                    early_exaggeration=st.session_state.tsne_early_exaggeration,
                    learning_rate=st.session_state.tsne_learning_rate,
                    random_state=st.session_state.tsne_random_state,
                )
        else:
            raise ValueError(
                f"visualization_type '{st.session_state.visualization_type}' not recognized"
            )
        st.session_state.dimension_reduction_complete = True


def run():
    if "dimension_reduction_complete" not in st.session_state:
        st.session_state.dimension_reduction_complete = False

    st.markdown("# Dimensional Reduction: UMAP/t-SNE")

    display_sidebar()

    if st.session_state.dimension_reduction_complete:
        if st.session_state.visualization_type == "UMAP":
            display_umap_results(st.session_state.normalized_adata)
        elif st.session_state.visualization_type == "t-SNE":
            display_tsne_results(st.session_state.normalized_adata)
        else:
            raise ValueError(
                f"visualization_type '{st.session_state.visualization_type}' not recognized"
            )


run()
