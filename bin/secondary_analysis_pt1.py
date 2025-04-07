#!/usr/bin/env python3

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path
from scipy.stats import median_abs_deviation

import anndata
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import scanpy as sc
import uuid

def add_patient_metadata(obs, uuids_df):
    merged = uuids_df.merge(obs, left_on="uuid", right_on="dataset", how="inner")
    merged = merged.set_index(obs.index)
    merged = merged.drop(columns=["Unnamed: 0"], errors='ignore')
    merged = merged.fillna(np.nan)
    merged["age"] = pd.to_numeric(merged["age"], errors='coerce')
    obs = obs.loc[:, ~obs.columns.str.contains("^Unnamed")]
    return merged


def add_file_sizes(data_product_metadata, raw_size):
    data_product_metadata["Raw File Size"] = raw_size


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def main(
    raw_h5ad_file: Path,
    uuids_tsv: Path,
    data_product_metadata: Path,
    tissue: str = None,
):
    raw_output_file_name = f"{tissue}_raw.h5ad" if tissue else "rna_raw.h5ad"
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )

    adata = anndata.read_h5ad(raw_h5ad_file)
    dataset_info = pd.read_csv(uuids_tsv, sep="\t")
    # annotated_obs = add_patient_metadata(adata.obs, dataset_info)
    # adata.obs = annotated_obs
    print("Writing raw data product")
    print(adata.obs_keys())
    adata.write(raw_output_file_name)

    raw_file_size = os.path.getsize(raw_output_file_name)
    with open(data_product_metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Data Product UUID"]
    add_file_sizes(metadata, raw_file_size)
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)

    print("Processing data product")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    # try new cell filtering method
    # mitochondrial genes
    adata.var["mt"] = adata.var["hugo_symbol"].str.startswith("MT-")
    adata.var["mt"] = adata.var["mt"].fillna(False)
    # ribosomal genes
    adata.var["ribo"] = adata.var["hugo_symbol"].str.startswith(("RPS", "RPL"))
    adata.var["ribo"] = adata.var["ribo"].fillna(False)
    # hemoglobin genes.
    adata.var["hb"] = adata.var["hugo_symbol"].str.contains("^HB[^(P)]")
    adata.var["hb"] = adata.var["hb"].fillna(False)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")  

    sc.pp.filter_genes(adata, min_cells=3)

    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["unscaled"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=50)
    sc.tl.embedding_density(adata, basis="pca")  

    adata.write(processed_output_file_name)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("uuids_file")
    p.add_argument("data_product_metadata")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.uuids_file, args.data_product_metadata, args.tissue)
