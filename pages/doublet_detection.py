import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.plotting import display_doublet_detection_info
from utils.preprocessing import detect_doublets
from utils.utils import set_downstream_pages_to_not_complete


@dataclass
class DoubletPageState:
    # user selections when "Detect Doublets" button was last clicked
    doublet_algorithm: Union[str, None] = None
    scrublet_seed: Union[int, None] = None
    scrublet_n_prin_comps: Union[int, None] = None
    scrublet_expected_doublet_rate: Union[float, None] = None

    # current user selections
    user_sel_doublet_algorithm: Union[str, None] = None
    user_sel_scrublet_seed: Union[int, None] = None
    user_sel_scrublet_n_prin_comps: Union[int, None] = None
    user_sel_scrublet_expected_doublet_rate: Union[float, None] = None

    # parameters not selected by user
    filtered_adata: Union[AnnData, None] = None
    doublet_detection_complete: bool = False
    doublet_step_complete: bool = False

    def reset(self):
        self.doublet_algorithm = None
        self.scrublet_seed = None
        self.scrublet_n_prin_comps = None
        self.scrublet_expected_doublet_rate = None

    def update(self):
        self.doublet_algorithm = self.user_sel_doublet_algorithm
        if self.doublet_algorithm == "Scrublet":
            self.scrublet_seed = self.user_sel_scrublet_seed
            self.scrublet_n_prin_comps = self.user_sel_scrublet_n_prin_comps
            self.scrublet_expected_doublet_rate = self.user_sel_scrublet_expected_doublet_rate
        else:
            raise ValueError(f"doublet_algorithm '{self.doublet_algorithm}' not recognized")


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            try:
                self.state = DoubletPageState(
                    filtered_adata=st.session_state.quality_control.filtered_adata
                )
            except Exception:
                self.state = DoubletPageState()
        self.detect_doublets_clicked = False
        self.remove_doublets_clicked = False

    def display_algorithm_advanced_options(self, algorithm):
        if algorithm == "Scrublet":
            with st.sidebar.expander("Advanced Options", expanded=False):
                self.state.user_sel_scrublet_seed = st.number_input(
                    "Random seed",
                    min_value=0,
                    step=1,
                    value=self.state.scrublet_seed if self.state.scrublet_seed else 0,
                )
                self.state.user_sel_scrublet_n_prin_comps = st.number_input(
                    "Num Principal Components",
                    min_value=1,
                    value=(
                        self.state.scrublet_n_prin_comps if self.state.scrublet_n_prin_comps else 30
                    ),
                )
                self.state.user_sel_scrublet_expected_doublet_rate = st.number_input(
                    "Expected Doublet Rate",
                    min_value=0.0,
                    max_value=1.0,
                    value=(
                        self.state.scrublet_expected_doublet_rate
                        if self.state.scrublet_expected_doublet_rate
                        else 0.05
                    ),
                    step=0.01,
                )
        elif algorithm == "Vaeda":
            raise NotImplementedError("Vaeda algorithm not implemented")
        else:
            raise ValueError(f"algorithm '{algorithm}' not recgonized")

    def display_sidebar(self):
        # Algorithm selection and advanced parameters
        algorithm_options = ["Scrublet"]
        self.state.user_sel_doublet_algorithm = st.sidebar.selectbox(
            "Doublet Detection Algorithm",
            options=algorithm_options,
            index=(
                algorithm_options.index(self.state.doublet_algorithm)
                if self.state.doublet_algorithm
                else 0
            ),
            help="The algorithm to use for detecting doublets",
        )
        self.display_algorithm_advanced_options(self.state.user_sel_doublet_algorithm)

        # Detect Doublets button
        self.detect_doublets_clicked = st.sidebar.button(
            "Detect Doublets",
            help="Run detect doublets",
        )

        # Remove Doublets button
        self.remove_doublets_clicked = st.sidebar.button(
            "Remove Doublets",
            help="Remove detected doublets",
            disabled=not self.state.doublet_detection_complete,
        )

    def remove_doublets(self, adata, algorithm="Scrublet"):
        if algorithm == "Scrublet":
            self.state.doublets_removed_adata = adata[~adata.obs["predicted_doublet"]].copy()
        elif algorithm == "Vaeda":
            raise NotImplementedError("Vaeda algorithm not implemented")
        else:
            raise ValueError("algorithm '{algorithm}' not recgonized")

        # save state to global st.session_state
        self.state.doublet_step_complete = True
        self.save_to_session_state()

    def run_detect_doublets(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        with st.spinner("Doublet detection running.  Please allow 2-3 minutes."):
            detect_doublets(
                self.state.filtered_adata,
                self.state,
                self.state.doublet_algorithm,
            )

        # save state to global st.session_state
        self.state.doublet_detection_complete = True
        self.state.doublet_step_complete = False
        self.save_to_session_state()

    def display_doublet_info(self):
        if self.state.doublet_detection_complete:
            st.markdown("## Before doublet removal")
            display_doublet_detection_info(
                self.state.filtered_adata, algorithm=self.state.doublet_algorithm
            )

            st.markdown("## After doublet removal")
            if self.state.doublet_step_complete:
                display_doublet_detection_info(
                    self.state.doublets_removed_adata,
                    algorithm=self.state.doublet_algorithm,
                )
            else:
                st.markdown(
                    "⚠️ No doublets have been removed.  To remove doublets, click the *Remove "
                    "Doublets* button."
                )
        else:
            st.markdown(
                "⚠️ Doublet detection has not been run.  Select the desired parameters on the"
                " left and click *Detect Doublets* to continue."
            )

    def save_to_session_state(self):
        st.session_state.doublet_detection = self.state

    def run(self):
        text = """# Doublet Detection

The Doublet Detection step identifies likely doublets and removes them from downstream analysis.

Doublets are two cells that were captured in the same GEM and were sequenced using the same cellular barcode.  Uncorrected, the presence of doublets in the dataset can be problematic in downstream steps.  For example, two cells with different expression profiles that are combined under the same cellular barcode may incorrectly appear as a unique subpopulation of cells or as part of a transition state. 
"""
        st.markdown(text)
        page_step_number = st.session_state.page_completion_order.index("doublet_detection")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("⚠️ Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()

        # "Detect Doublets" button clicked
        if self.detect_doublets_clicked:
            self.run_detect_doublets()

            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number - 1
            set_downstream_pages_to_not_complete("doublet_detection")

            # this allows the page to properly enable the "Remove Doublets" button
            st.rerun()

        # "Remove Doublets" button clicked
        if self.remove_doublets_clicked:
            self.remove_doublets(self.state.filtered_adata, algorithm=self.state.doublet_algorithm)
            st.write("Doublets have been removed!")

            # update with furthest step completed
            st.session_state.furthest_step_number_completed = page_step_number

        # displaying the sidebar comes after button click handling so that it properly enables the
        # the "Remove Doublets" button after "Detect Doublets" has been run

        # Display doublet information
        self.display_doublet_info()


if "doublet_detection" in st.session_state:
    page = Page(page_state=st.session_state.doublet_detection)
else:
    page = Page()
page.run()
