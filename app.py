from pathlib import Path

import streamlit as st

# using "wide" allows more information to be displayed, but st.pyplot()
# fills the entire container width, which can be a little unwieldy
# st.set_page_config(layout="centered")
st.set_page_config(layout="wide")

dir_path = Path(__file__).parent


def main():
    st.session_state.page_completion_order = [
        "introduction",
        "load_data",
        "quality_control",
        "doublet_detection",
        "normalization",
        "feature_selection",
        "pca",
        "projection",
        "clustering",
    ]
    if "furthest_step_number_completed" not in st.session_state:
        st.session_state.furthest_step_number_completed = 0
    page = st.navigation(
        [
            st.Page(f"{dir_path}/pages/introduction.py", title="Introduction"),
            st.Page(f"{dir_path}/pages/load_data.py", title="Load Data"),
            st.Page(f"{dir_path}/pages/quality_control.py", title="Quality Control"),
            st.Page(f"{dir_path}/pages/doublet_detection.py", title="Doublet Detection"),
            st.Page(f"{dir_path}/pages/normalization.py", title="Normalization"),
            st.Page(f"{dir_path}/pages/feature_selection.py", title="Feature Selection"),
            st.Page(f"{dir_path}/pages/dim_reduction_pca.py", title="Dim Reduction: PCA"),
            st.Page(
                f"{dir_path}/pages/dim_reduction_umap_tsne.py", title="Dim Reduction: UMAP / t-SNE"
            ),
            st.Page(f"{dir_path}/pages/clustering.py", title="Clustering"),
        ]
    )
    page.run()


if __name__ == "__main__":
    main()
