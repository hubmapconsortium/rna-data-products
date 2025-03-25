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


def add_cell_counts(data_product_metadata, cell_counts, total_cell_count):
    with open(data_product_metadata, "r") as json_file:
        metadata = json.load(json_file)
    metadata["Processed Cell Type Counts"] = cell_counts
    metadata["Processed Total Cell Count"] = total_cell_count
    return metadata


def add_file_sizes(data_product_metadata, processed_size):
    data_product_metadata["Processed File Size"] = processed_size
    uuid = data_product_metadata["Data Product UUID"]
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(data_product_metadata, outfile)


def main(h5ad_file: Path, data_product_metadata: Path, tissue: str=None):
    adata = anndata.read_h5ad(h5ad_file)
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )
    total_cell_count = adata.obs.shape[0]
    with open(data_product_metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Data Product UUID"]
    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
    sc.tl.umap(adata)

    # leiden clustering
    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, "leiden")

    if "predicted_label" in adata.obs_keys():

        non_na_values = adata.obs.predicted_label.dropna()
        counts_dict = non_na_values.value_counts().to_dict()
        keep_cell_types = [
            cell_type for cell_type in counts_dict if counts_dict[cell_type] > 1
        ]
        adata_filter = adata[adata.obs.predicted_label.isin(keep_cell_types)]
        sc.tl.rank_genes_groups(
            adata_filter, "predicted_label", key_added="rank_genes_groups_cell_types"
        )
        adata.uns = adata_filter.uns

    if "predicted_label" in adata.obs_keys():
        cell_type_counts = adata.obs["predicted_label"].value_counts().to_dict()
        adata.uns["cell_type_counts"] = cell_type_counts
        metadata = add_cell_counts(
            data_product_metadata, cell_type_counts, total_cell_count
        )
    else:
        cell_type_counts = {}
        metadata = add_cell_counts(
            data_product_metadata, cell_type_counts, total_cell_count
        )

    with plt.rc_context():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig(f"{uuid}.png" if tissue else "rna.png")

    print(f"Writing {processed_output_file_name}")
    adata.write(processed_output_file_name)
    processed_file_size = os.path.getsize(processed_output_file_name)
    add_file_sizes(metadata, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("data_product_metadata", type=Path)
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.data_product_metadata, args.tissue)