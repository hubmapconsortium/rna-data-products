#!/usr/bin/env python3

from argparse import ArgumentParser
from datetime import datetime
from pathlib import Path

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


def main(
    raw_h5ad_file: Path,
    tissue: str = None,
):
    raw_output_file_name = f"{tissue}_raw.h5ad" if tissue else "rna_raw.h5ad"
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )

    adata = anndata.read_h5ad(raw_h5ad_file)
    print("Processing data product")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
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
    # p.add_argument("uuids_file")
    # p.add_argument("data_product_metadata")
    # p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.tissue)
