import scanpy as sc
import streamlit as st


def run_quality_control(adata):
    # annotate genes: mitochondrial ("mt"), ribosomal ("rb"), hemoglobin ("hb")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    # run QC
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    st.session_state.quality_control_complete = True
