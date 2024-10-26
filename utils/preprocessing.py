import scanpy as sc
import streamlit as st


def run_quality_control(adata):
    # annotate genes: mitochondrial ("mt"), ribosomal ("rb"), hemoglobin ("hb")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    # run QC
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)


def detect_doublets(adata, session_state, algorithm="Scrublet"):
    if algorithm == "Scrublet":
        seed = session_state.scrublet_seed
        n_prin_comps = session_state.scrublet_n_prin_comps
        expected_doublet_rate = session_state.scrublet_expected_doublet_rate
        detect_doublets_scrublet(
            adata, seed=seed, n_prin_comps=n_prin_comps, expected_doublet_rate=expected_doublet_rate
        )
    elif algorithm == "Vaeda":
        raise NotImplementedError("Vaeda algorithm not implemented")
    else:
        raise ValueError(f"algorithm '{algorithm}' not recgonized")


def detect_doublets_scrublet(adata, seed=0, n_prin_comps=30, expected_doublet_rate=0.05):
    sc.pp.scrublet(
        adata,
        batch_key="sample",
        random_state=seed,
        n_prin_comps=n_prin_comps,
        expected_doublet_rate=expected_doublet_rate,
    )
