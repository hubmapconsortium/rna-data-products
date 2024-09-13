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


def annotate_h5ad(obs, uuids_df: pd.DataFrame):
    uuids_df["uuid"] = uuids_df["uuid"].astype(str)
    obs["dataset"] = obs["dataset"].astype(str)
    uuids_dict = uuids_df.set_index('uuid').to_dict(orient='index')
    for index, row in obs.iterrows():
        dataset_value = row['dataset']
        if dataset_value in uuids_dict:
            for key, value in uuids_dict[dataset_value].items():
                obs.at[index, key] = value
    obs.drop(columns=['Unnamed: 0'])
    return obs


def add_cell_counts(data_product_metadata, cell_counts, total_cell_count):
    with open(data_product_metadata, 'r') as json_file:
        metadata = json.load(json_file)
    metadata["Processed Cell Type Counts"] = cell_counts
    metadata["Processed Total Cell Count"] = total_cell_count
    return metadata

def add_file_sizes(data_product_metadata, raw_size, processed_size):
    data_product_metadata["Raw File Size"] = raw_size
    data_product_metadata["Processed File Size"] = processed_size
    uuid = data_product_metadata["Data Product UUID"]
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(data_product_metadata, outfile)


def main(raw_h5ad_file: Path, data_product_metadata: Path, uuids_tsv: Path, tissue: str = None):
    raw_output_file_name = f"{tissue}_raw" if tissue else "rna_raw"
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )
    
    raw_file_size = os.path.getsize(raw_h5ad_file)

    adata = anndata.read_h5ad(raw_h5ad_file)
    dataset_info = pd.read_csv(uuids_tsv, sep="\t")
    annotated_obs = annotate_h5ad(adata.obs, dataset_info)
    adata.obs = annotated_obs
    print("Writing raw data product")
    print(adata.obs_keys())
    adata.write_h5ad(f"{raw_output_file_name}.h5ad")

    uuid = str(adata.uns["uuid"])
    
    print("Processing data product")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    total_cell_count = adata.obs.shape[0]
    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["unscaled"] = adata.X.copy()
    #    sc.pp.combat(adata, "dataset")
    sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)

    sc.tl.umap(adata)

    # leiden clustering
    sc.tl.leiden(adata)

    if "predicted_label" in adata.obs_keys():

        non_na_values = adata.obs.predicted_label.dropna()
        counts_dict = non_na_values.value_counts().to_dict()
        keep_cell_types = [
            cell_type for cell_type in counts_dict if counts_dict[cell_type] > 1
        ]
        adata_filter = adata[adata.obs.predicted_label.isin(keep_cell_types)]
        # Filter out cell types with only one cell for this analysis
        # sc.pp.filter_genes(adata_filter, min_cells=3)
        sc.tl.rank_genes_groups(adata_filter, "predicted_label")
        adata.uns = adata_filter.uns
    
    if "predicted_label" in adata.obs_keys():
        cell_type_counts = (adata.obs["predicted_label"].value_counts().to_dict())
        adata.uns["cell_type_counts"] = cell_type_counts
        metadata = add_cell_counts(data_product_metadata, cell_type_counts, total_cell_count)
    else:
        cell_type_counts = {}
        metadata = add_cell_counts(data_product_metadata, cell_type_counts, total_cell_count)

    with plt.rc_context():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig(f"{uuid}.png" if tissue else "rna.png")

    print(f"Writing {processed_output_file_name}")
    adata.write(processed_output_file_name)
    processed_file_size = os.path.getsize(processed_output_file_name)
    add_file_sizes(metadata, raw_file_size, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("uuids_file")
    p.add_argument("data_product_metadata", type=Path)
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.data_product_metadata, args.uuids_file, args.tissue)
