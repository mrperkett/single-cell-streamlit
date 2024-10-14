import numpy as np
import pandas as pd


def get_qc_stats_dataframe(adata):
    # total UMIs for each gene
    total_umis_by_gene = adata.var["total_counts"].astype(int)
    idx = total_umis_by_gene > 0
    min_umis_per_gene = np.min(total_umis_by_gene[idx])
    median_umis_per_gene = np.median(total_umis_by_gene[idx])
    mean_umis_per_gene = np.mean(total_umis_by_gene[idx])
    max_umis_per_gene = np.max(total_umis_by_gene[idx])
    stdev_umis_per_gene = np.std(total_umis_by_gene[idx])

    # total UMIs for each cell
    total_umis_by_cell = adata.obs["total_counts"].astype(int)
    idx = total_umis_by_cell > 0
    min_umis_per_cell = np.min(total_umis_by_cell[idx])
    median_umis_per_cell = np.median(total_umis_by_cell[idx])
    mean_umis_per_cell = np.mean(total_umis_by_cell[idx])
    max_umis_per_cell = np.max(total_umis_by_cell[idx])
    stdev_umis_per_cell = np.std(total_umis_by_cell[idx])

    # total genes for each cell
    total_genes_by_cell = adata.obs["n_genes_by_counts"]
    idx = total_genes_by_cell > 0
    min_genes_per_cell = np.min(total_genes_by_cell[idx])
    median_genes_per_cell = np.median(total_genes_by_cell[idx])
    mean_genes_per_cell = np.mean(total_genes_by_cell[idx])
    max_genes_per_cell = np.max(total_genes_by_cell[idx])
    stdev_genes_per_cell = np.std(total_genes_by_cell[idx])

    # totals over
    total_umi_count = adata.X.astype(int).sum()
    total_genes_identified = np.count_nonzero(total_umis_by_gene)

    data = [
        ["Total UMI Count", total_umi_count],
        ["Total Genes Identified", total_genes_identified],
        ["Median UMIs per Gene", median_umis_per_gene],
        ["Mean UMIs per Gene", mean_umis_per_gene],
        ["Median UMIs per Cell", median_umis_per_cell],
        ["Mean UMIs per Cell", mean_umis_per_cell],
    ]
    df_general = pd.DataFrame(data, columns=["Label", "Value"])

    data = [
        ["Min UMIs per Gene", min_umis_per_gene],
        ["Median UMIs per Gene", median_umis_per_gene],
        ["Mean UMIs per Gene", mean_umis_per_gene],
        ["Max UMIs per Gene", max_umis_per_gene],
        ["Stev UMIs per Gene", stdev_umis_per_gene],
    ]
    df_umis_per_gene = pd.DataFrame(data, columns=["Label", "Value"])

    data = [
        ["Min UMIs per Cell", min_umis_per_cell],
        ["Median UMIs per Cell", median_umis_per_cell],
        ["Mean UMIs per Cell", mean_umis_per_cell],
        ["Max UMIs per Cell", max_umis_per_cell],
        ["Stdev UMIs per Cell", stdev_umis_per_cell],
    ]
    df_umis_per_cell = pd.DataFrame(data, columns=["Label", "Value"])

    data = [
        ["Min genes per Cell", min_genes_per_cell],
        ["Median genes per Cell", median_genes_per_cell],
        ["Mean genes per Cell", mean_genes_per_cell],
        ["Max genes per Cell", max_genes_per_cell],
        ["Stdev genes per Cell", stdev_genes_per_cell],
    ]
    df_genes_per_cell = pd.DataFrame(data, columns=["Label", "Value"])

    return df_general, df_umis_per_gene, df_umis_per_cell, df_genes_per_cell
