import copy
from dataclasses import dataclass
from typing import Union

import streamlit as st
from anndata import AnnData

from utils.analysis import filter_adata
from utils.plotting import display_qc_info_comparison
from utils.utils import set_downstream_pages_to_not_complete


@dataclass
class QualityControlPageState:
    # user selections when "Run Dimensional Reduction" button was last clicked
    min_allowed_genes_in_cell: Union[int, None] = None
    min_allowed_cells_with_gene: Union[int, None] = None
    max_allowed_percent_mt: Union[float, None] = None

    # current user selections
    user_sel_min_allowed_genes_in_cell: Union[int, None] = None
    user_sel_min_allowed_cells_with_gene: Union[int, None] = None
    user_sel_max_allowed_percent_mt: Union[float, None] = None

    # parameters not selected by user
    adata: Union[AnnData, None] = None
    filtered_adata: Union[AnnData, None] = None
    quality_control_complete: bool = False

    def reset(self):
        self.min_allowed_genes_in_cell = None
        self.min_allowed_cells_with_gene = None
        self.max_allowed_percent_mt = None

    def update(self):
        self.min_allowed_genes_in_cell = self.user_sel_min_allowed_genes_in_cell
        self.min_allowed_cells_with_gene = self.user_sel_min_allowed_cells_with_gene
        self.max_allowed_percent_mt = self.user_sel_max_allowed_percent_mt


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            try:
                self.state = QualityControlPageState(adata=st.session_state.load_data.adata)
            except Exception:
                self.state = QualityControlPageState()

    def run_filter(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()

        self.state.filtered_adata = filter_adata(
            self.state.adata,
            self.state.min_allowed_genes_in_cell,
            self.state.min_allowed_cells_with_gene,
            self.state.max_allowed_percent_mt,
        )

        # save state to global st.session_state
        self.state.quality_control_complete = True
        self.save_to_session_state()

    def display_sidebar(self):
        self.state.user_sel_min_allowed_genes_in_cell = st.sidebar.number_input(
            "Min allowed genes in a cell",
            min_value=1,
            value=(
                self.state.min_allowed_genes_in_cell
                if self.state.min_allowed_genes_in_cell
                else "min"
            ),
            step=1,
            help="Cells that have less than this number of genes will be removed",
        )

        self.state.user_sel_min_allowed_cells_with_gene = st.sidebar.number_input(
            "Min allowed cells with gene",
            min_value=1,
            value=(
                self.state.min_allowed_cells_with_gene
                if self.state.min_allowed_cells_with_gene
                else "min"
            ),
            step=1,
            help="Genes that are identified in less than this number of cells will be removed",
        )

        self.state.user_sel_max_allowed_percent_mt = st.sidebar.slider(
            "Max allowed % Mitochondrial",
            min_value=0.0,
            max_value=100.0,
            value=self.state.max_allowed_percent_mt if self.state.max_allowed_percent_mt else 100.0,
            step=1.0,
            help="Cells where the percent Mitochondrial UMIs is greater than this number will be"
            "removed",
        )

        self.apply_filters_button_clicked = st.sidebar.button(
            "Apply Filters",
            help="Apply quality control filters to data",
        )

    def save_to_session_state(self):
        st.session_state.quality_control = self.state

    def run(self):
        text = """# Quality Control

The Quality Control step filters out undesirable data based on user-defined thresholds.  This includes

- barcodes that have fewer than *min allowed genes in a cell*
- barcodes with more than *max allowed % mitochondrial* mitochondrial UMIs (i.e. unique identified transcripts)
- genes that occur in fewer than *min allowed cells with gene* barcodes
        
Each of these cutoffs is meant to filter out low-quality cells from the dataset.  Low-quality barcodes may result from dying cells, which have ruptured membranes and have released cytoplasmic mRNA leaving only mitochondrial mRNA.

Care should be taken when defining cutoffs since it can be difficult to distinguish between outliers due to meaningful biological differences and a systematic error in the experiment.  It is generally recommended to be as permissive as possible to avoid filtering out outlier cell sub-populations.
"""
        st.markdown(text)
        page_step_number = st.session_state.page_completion_order.index("quality_control")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("⚠️ Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.display_sidebar()
        if self.apply_filters_button_clicked:
            self.run_filter()

            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number
            set_downstream_pages_to_not_complete("quality_control")
        display_qc_info_comparison(self.state.adata, self.state.filtered_adata)


if "quality_control" in st.session_state:
    page = Page(page_state=st.session_state.quality_control)
else:
    page = Page()
page.run()
