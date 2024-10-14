import streamlit as st

from utils.plotting import display_doublet_detection_info
from utils.preprocessing import detect_doublets


def display_algorithm_options(algorithm):
    if algorithm == "Scrublet":
        with st.sidebar.expander("Advanced Options", expanded=False):
            st.session_state.scrublet_seed = st.number_input(
                "Random seed", min_value=0, step=1, value=0
            )
            st.session_state.scrublet_n_prin_comps = st.number_input(
                "Num Principal Components", min_value=1, value=30
            )
            st.session_state.scrublet_expected_doublet_rate = st.number_input(
                "Expected Doublet Rate", min_value=0.0, max_value=1.0, value=0.05, step=0.01
            )
    elif algorithm == "Vaeda":
        raise NotImplementedError("Vaeda algorithm not implemented")
    else:
        raise ValueError("algorithm '{algorithm}' not recgonized")


def display_sidebar():
    # Algorithm selection and advanced parameters
    algorithm_options = ["Scrublet"]
    st.session_state.doublet_algorithm = st.sidebar.selectbox(
        "Doublet Detection Algorithm",
        options=algorithm_options,
        help="The algorithm to use for detecting doublets",
    )

    display_algorithm_options(st.session_state.doublet_algorithm)

    # Detect Doublets button
    st.session_state.detect_doublets_clicked = st.sidebar.button(
        "Detect Doublets",
        help="Run detect doublets",
    )

    # Remove Doublets button
    st.session_state.remove_doublets_clicked = st.sidebar.button(
        "Remove Doublets",
        help="Remove detected doublets",
        disabled=not st.session_state.doublets_detection_complete,
    )


def remove_doublets(adata, algorithm="Scrublet"):
    if algorithm == "Scrublet":
        doublets_removed_adata = adata[~adata.obs["predicted_doublet"]].copy()
    elif algorithm == "Vaeda":
        raise NotImplementedError("Vaeda algorithm not implemented")
    else:
        raise ValueError("algorithm '{algorithm}' not recgonized")
    return doublets_removed_adata


def run():

    if "doublets_detection_complete" not in st.session_state:
        st.session_state.doublets_detection_complete = False

    display_sidebar()

    # action when "Detect Doublets" is clicked
    if st.session_state.detect_doublets_clicked:
        with st.spinner("Doublet detection running.  Please allow 2-3 minutes."):
            detect_doublets(
                st.session_state.filtered_adata,
                st.session_state,
                st.session_state.doublet_algorithm,
            )
        st.session_state.doublets_detection_complete = True
        st.session_state.removed_doublets_clicked = False
        st.session_state.doublets_removed_adata = None

    # action when "Remove Doublets" is clicked
    if st.session_state.remove_doublets_clicked:
        st.session_state.doublets_removed_adata = remove_doublets(
            st.session_state.filtered_adata, algorithm=st.session_state.doublet_algorithm
        )
        st.write("Doublets have been removed!")

    # Display doublet information
    st.markdown("# Doublet Detection")
    if st.session_state.doublets_detection_complete:

        st.markdown("## Before doublet removal")
        display_doublet_detection_info(
            st.session_state.filtered_adata, algorithm=st.session_state.doublet_algorithm
        )

        st.markdown("## After doublet removal")
        if "doublets_removed_adata" in st.session_state and st.session_state.doublets_removed_adata:
            display_doublet_detection_info(
                st.session_state.doublets_removed_adata,
                algorithm=st.session_state.doublet_algorithm,
            )
        else:
            st.markdown(
                "No doublets have been removed.  To remove doublets, click the *Remove Doublets* button."
            )


run()
