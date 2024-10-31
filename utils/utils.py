import streamlit as st


def get_downstream_pages(page_name):
    if page_name not in st.session_state.page_completion_order:
        raise ValueError(f"page_name '{page_name}' not recognized")

    page_step_number = st.session_state.page_completion_order.index(page_name)
    downstream_page_names = [
        page_name for page_name in st.session_state.page_completion_order[page_step_number + 1 :]
    ]

    return downstream_page_names


def set_page_to_not_complete(page_name):
    # if the page has not been saved to the st.session_state, then there is nothing to do
    if page_name not in st.session_state:
        return

    # otherwise, each page has a unique way of being reset
    if page_name == "quality_control":
        st.session_state.quality_control.quality_control_complete = False
        st.session_state.quality_control.filtered_adata = None
    elif page_name == "doublet_detection":
        st.session_state.doublet_detection.doublet_detection_complete = False
        st.session_state.doublet_detection.doublet_step_complete = False
    elif page_name == "normalization":
        st.session_state.normalization.run_normalization_complete = False
    elif page_name == "feature_selection":
        st.session_state.feature_selection.feature_selection_complete = False
    elif page_name == "pca":
        st.session_state.pca.run_pca_complete = False
    elif page_name == "projection":
        st.session_state.projection.projection_complete = False
    elif page_name == "clustering":
        st.session_state.clustering.clustering_complete = False
    else:
        raise ValueError(f"page_name '{page_name}' not recognized")


def set_downstream_pages_to_not_complete(page_name):
    downstream_page_names = get_downstream_pages(page_name)

    for page_name in downstream_page_names:
        set_page_to_not_complete(page_name)
