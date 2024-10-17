import scanpy as sc


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


def run_pca(adata, n_comps, use_highly_variable, random_state):
    sc.tl.pca(
        adata, n_comps=n_comps, use_highly_variable=use_highly_variable, random_state=random_state
    )


def run_umap(adata, min_dist, spread, random_state, alpha):
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=random_state, alpha=alpha)


def run_tsne(adata, n_pcs, perplexity, early_exaggeration, learning_rate, random_state):
    sc.tl.tsne(
        adata,
        n_pcs=n_pcs,
        perplexity=perplexity,
        early_exaggeration=early_exaggeration,
        learning_rate=learning_rate,
        random_state=random_state,
    )


def run_clustering(adata, page_state):
    method = page_state.clustering_method
    if method == "Leiden":
        sc.tl.leiden(
            adata,
            flavor="igraph",
            resolution=page_state.leiden_resolution,
            n_iterations=page_state.leiden_n_iterations,
            random_state=page_state.leiden_random_state,
        )
    elif method == "Louvain":
        sc.tl.louvain(
            adata,
            resolution=page_state.louvain_resolution,
            random_state=page_state.louvain_random_state,
        )
    else:
        raise ValueError(f"method '{method}' not recognized")
