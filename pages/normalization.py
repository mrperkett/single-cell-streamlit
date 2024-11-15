import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import run_normalization
from utils.plotting import display_get_umi_distribution_plot
from utils.utils import set_downstream_pages_to_not_complete


@dataclass
class NormalizationPageState:
    # user selections when "Run Dimensional Reduction" button was last clicked
    exclude_highly_expressed: Union[bool, None] = None
    max_fraction: Union[float, None] = None
    target_sum_category: Union[str, None] = None
    target_sum: Union[int, None] = None

    # current user selections
    user_sel_exclude_highly_expressed: Union[bool, None] = None
    user_sel_max_fraction: Union[float, None] = None
    user_sel_target_sum_category: Union[str, None] = None
    user_sel_target_sum: Union[int, None] = None

    # parameters not selected by user
    run_normalization_complete: bool = False
    doublets_removed_adata: Union[AnnData, None] = None

    def reset(self):
        self.exclude_highly_expressed = None
        self.max_fraction = None
        self.target_sum_category = None
        self.target_sum = None
        self.run_normalization_complete = False

    def update(self):
        self.exclude_highly_expressed = self.user_sel_exclude_highly_expressed
        self.max_fraction = self.user_sel_max_fraction
        self.target_sum_category = self.user_sel_target_sum_category
        self.target_sum = self.user_sel_target_sum


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            try:
                self.state = NormalizationPageState(
                    doublets_removed_adata=st.session_state.doublet_detection.doublets_removed_adata
                )
            except Exception:
                self.state = NormalizationPageState()

    def display_sidebar(self):
        self.state.user_sel_exclude_highly_expressed = st.sidebar.checkbox(
            "Exclude Highly Expressed Genes",
            value=(
                self.state.exclude_highly_expressed if self.state.exclude_highly_expressed else True
            ),
            help="From Scanpy documentation: "
            "Exclude (very) highly expressed genes for the computation of the normalization factor"
            " (size factor) for each cell. A gene is considered highly expressed, if it has more "
            "than max_fraction of the total counts in at least one cell. The not-excluded genes"
            "will sum up to target_sum. Providing this argument when adata.X is a Array will incur"
            " blocking .compute() calls on the array. [From scanpy documentation]",
        )
        self.state.user_sel_max_fraction = st.sidebar.number_input(
            "Max Fraction",
            min_value=0.0,
            max_value=1.0,
            value=self.state.max_fraction if self.state.max_fraction else 0.05,
            step=0.01,
            help="From Scanpy documentation: "
            "If exclude_highly_expressed=True, consider cells as highly expressed that have more "
            "counts than max_fraction of the original total counts in at least one cell.",
            disabled=not self.state.exclude_highly_expressed,
        )
        options_list = ["Median", "Custom"]
        self.state.user_sel_target_sum_category = st.sidebar.selectbox(
            "Target Sum",
            options=options_list,
            index=(
                options_list.index(self.state.target_sum_category)
                if self.state.target_sum_category
                else 0
            ),
            help="From Scanpy documentation: "
            "If None, after normalization, each observation (cell) has a total count equal to the"
            " median of total counts for observations (cells) before normalization.",
        )

        self.state.user_sel_target_sum = st.sidebar.number_input(
            "Target Sum Custom Value",
            value=self.state.target_sum if self.state.target_sum else 10000,
            disabled=self.state.user_sel_target_sum_category != "Custom",
        )

        self.run_normalization_button_clicked = st.sidebar.button(
            "Run Normalization",
        )

    def run_normalization(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        if self.state.target_sum_category == "Median":
            target_sum = None
        else:
            target_sum = self.state.target_sum
        self.state.normalized_adata = run_normalization(
            self.state.doublets_removed_adata,
            self.state.exclude_highly_expressed,
            self.state.max_fraction,
            target_sum,
        )

        # save state to global st.session_state
        self.state.run_normalization_complete = True
        self.save_to_session_state()

    def save_to_session_state(self):
        st.session_state.normalization = self.state

    def display_normalization(self):

        if self.state.run_normalization_complete:
            st.write("Data has been normalized!")
            col1, col2 = st.columns(2)

            with col1:
                st.write("## Unnormalized")
                display_get_umi_distribution_plot(self.state.doublets_removed_adata)
            with col2:
                st.write("## Normalized")
                display_get_umi_distribution_plot(self.state.normalized_adata)
        else:
            st.markdown(
                "⚠️ Normalization has not been run.  Select the desired parameters on the"
                " left and click *Run Normalization*"
            )

    def run(self):
        text = """# Normalization

The Normalization step transforms the data to remove variable sample effects.

This is important to do since the data may vary purely due to various in the experimental procedure (e.g. transcript capture, reverse transcription, sequencing, etc).  This step first scales the total number of UMIs for a cell to a user-specified value (default = median) and then carries out a *log(count+1)* transformation.
"""
        st.markdown(text)
        page_step_number = st.session_state.page_completion_order.index("normalization")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("⚠️ Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()
        if self.run_normalization_button_clicked:
            self.run_normalization()
            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number
            set_downstream_pages_to_not_complete("normalization")
        self.display_normalization()


if "normalization" in st.session_state:
    page = Page(page_state=st.session_state.normalization)
else:
    page = Page()
page.run()
