import copy
from dataclasses import dataclass
from typing import Union

import anndata as ad
import scanpy as sc
import streamlit as st
import yaml

from utils.plotting import display_qc_info
from utils.preprocessing import run_quality_control
from utils.utils import set_downstream_pages_to_not_complete


@dataclass
class LoadDataPageState:
    # user selections when "Load" button was last clicked
    dataset_name: Union[str, None] = None

    # current user selections
    user_sel_dataset_name: Union[str, None] = None

    # parameters not selected by user
    dataset_loaded: bool = False
    quality_control_complete: bool = False
    config: Union[str, None] = None
    dataset: Union[dict, None] = None

    def reset(self):
        self.dataset_name = None

    def update(self):
        self.dataset_name = self.user_sel_dataset_name


def load_config(config_file_path):
    with open(config_file_path, "r") as input_file:
        config = yaml.safe_load(input_file)
    return config


class Page:

    def __init__(self, page_state=None):
        if page_state:
            self.state = copy.copy(page_state)
        else:
            self.state = LoadDataPageState()

    def load_dataset(self):
        # Update state with current user selections
        self.state.reset()
        self.state.update()
        self.state.dataset = self.state.config["datasets"][self.state.dataset_name]

        num_samples = len(self.state.dataset["samples"])
        with st.status(f"Loading {num_samples} samples", expanded=True) as status:
            adatas = {}
            for sample_id, sample_path in self.state.dataset["samples"].items():
                st.write(f"Loading sample {sample_id}")
                sample_adata = sc.read_10x_h5(sample_path)
                sample_adata.var_names_make_unique()
                adatas[sample_id] = sample_adata

            st.write("Concatenating samples")
            adata = ad.concat(adatas, label="sample")
            adata.obs_names_make_unique()
        status.update(label="Load complete!", state="complete", expanded=False)

        # save state to global st.session_state
        self.state.adata = adata
        self.state.dataset_loaded = True
        self.state.quality_control_complete = False
        self.save_to_session_state()

    def display_sidebar(self):
        dataset_names = list(self.state.config["datasets"].keys())
        self.state.user_sel_dataset_name = st.sidebar.selectbox(
            "Dataset",
            options=dataset_names,
            index=(dataset_names.index(self.state.dataset_name) if self.state.dataset_name else 0),
            help="Select which dataset to load",
        )
        self.load_clicked = st.sidebar.button("Load")

    def save_to_session_state(self):
        st.session_state.load_data = self.state

    def run(self):
        st.markdown("# Load Data")
        page_step_number = st.session_state.page_completion_order.index("load_data")

        # If the previous step has not been completed, display a message to the user and return
        if st.session_state.furthest_step_number_completed < page_step_number - 1:
            st.write("⚠️ Please complete the previous step before running this step")
            return

        # Otherwise, run the page
        self.state.config = load_config("config.yaml")
        self.display_sidebar()
        if self.load_clicked:
            self.load_dataset()

            # update with furthest step completed and reset downstream pages to show not complete
            st.session_state.furthest_step_number_completed = page_step_number
            set_downstream_pages_to_not_complete("load_data")

        if self.state.dataset_loaded:
            st.write(f"**Loaded data**: {self.state.dataset_name}")
            description = f"*{self.state.dataset['description']}*"
            st.markdown(f"{description}")
            st.write("**Samples:**")
            samples_markdown_text = "\n".join(
                ["- " + sample for sample in self.state.dataset["samples"]]
            )
            st.markdown(samples_markdown_text)

        st.markdown("## Quality Control Plots")

        # Run Quality Control if the dataset has been loaded and it hasn't already been run
        if self.state.dataset_loaded and not self.state.quality_control_complete:
            with st.spinner("Running QC..."):
                run_quality_control(self.state.adata)
                self.state.quality_control_complete = True
                # save state to global st.session_state
                self.save_to_session_state()

        # Display Quality Control information after it has been run
        if self.state.quality_control_complete:
            display_qc_info(self.state.adata)
        else:
            st.markdown(
                "⚠️ Load Data has not yet been run. Select the desired dataset on the left"
                " and click *Load* to continue."
            )


if "load_data" in st.session_state:
    page = Page(page_state=st.session_state.load_data)
else:
    page = Page()
page.run()
