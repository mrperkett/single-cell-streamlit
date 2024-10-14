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
        raise ValueError("algorithm '{algorithm}' not recgonized")


def detect_doublets_scrublet(adata, seed=0, n_prin_comps=30, expected_doublet_rate=0.05):
    sc.pp.scrublet(
        adata,
        batch_key="sample",
        random_state=seed,
        n_prin_comps=n_prin_comps,
        expected_doublet_rate=expected_doublet_rate,
    )
    # import time
    # time.sleep(3)


def run_normalization(adata, exclude_highly_expressed, max_fraction, target_sum):
    normalized_adata = adata.copy()
    normalized_adata.layers["counts"] = normalized_adata.X.copy()
    sc.pp.normalize_total(
        normalized_adata,
        exclude_highly_expressed=exclude_highly_expressed,
        max_fraction=max_fraction,
        target_sum=target_sum,
    )

    # Log transform that data after adding one (i.e. log(counts + 1))
    sc.pp.log1p(normalized_adata)

    return normalized_adata


def mark_highly_variable_genes(adata, algorithm, **kwargs):
    if algorithm == "seurat":
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample", flavor=algorithm)
    elif algorithm == "cell_ranger":
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample", flavor=algorithm)
    elif algorithm == "seurat_v3":
        # Using layer="counts" since documentation states
        # Expects logarithmized data, except when flavor='seurat_v3'/'seurat_v3_paper', in which
        # count data is expected.
        sc.pp.highly_variable_genes(
            adata, n_top_genes=2000, batch_key="sample", flavor=algorithm, layer="counts"
        )
    else:
        raise ValueError(f"algorithm '{algorithm}' not recognized")
