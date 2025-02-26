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


def create_json(data_product_uuid, creation_time, raw_cell_count, processed_total_cell_count, raw_cell_type_counts, processed_cell_type_counts, raw_file_size, processed_file_size, tissue = None):
    bucket_url = f"https://hubmap-data-products.s3.amazonaws.com/{data_product_uuid}/"
    metadata = {
        "Data Product UUID": data_product_uuid,
        "Tissue": "Pan-Organ",
        "Assay": "rna",
        "Raw URL": bucket_url + f"{tissue}_raw.h5ad" if tissue else bucket_url + "RNA_raw.h5ad",
        "Processed URL": bucket_url + f"{tissue}_processed.h5ad" if tissue else bucket_url + "RNA_processed.h5ad",
        "Creation Time": creation_time,
        "Raw Total Cell Count": raw_cell_count,
        "Processed Total Cell Count": processed_total_cell_count,
        "Raw Cell Type Counts": raw_cell_type_counts,
        "Processed Cell Type Counts": processed_cell_type_counts,
        "Raw File Size": raw_file_size,
        "Processed File Size": processed_file_size

    }
    print("Writing metadata json")
    with open(f"{data_product_uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)


def main(
    raw_h5ad_file: Path,
    tissue: str = None,
):
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "RNA_processed.h5ad"
    )

    raw_file_size = os.path.getsize(raw_h5ad_file)

    adata = anndata.read_h5ad(raw_h5ad_file)

    data_product_uuid = str(uuid.uuid4())
    date_created = str(datetime.now())
    raw_cell_type_counts = adata.obs["predicted_label"].value_counts().to_dict()
    raw_total_cell_count = adata.obs.shape[0]

    print("Processing data product")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    processed_total_cell_count = adata.obs.shape[0]
    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["unscaled"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)

    sc.pp.pca(adata, n_comps=50)
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
    else:
        cell_type_counts = {}

    with plt.rc_context():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig(f"{data_product_uuid}.png")

    print(f"Writing {processed_output_file_name}")
    adata.write(processed_output_file_name)
    processed_file_size = os.path.getsize(processed_output_file_name)
    create_json(data_product_uuid, date_created, raw_total_cell_count, processed_total_cell_count, raw_cell_type_counts, cell_type_counts, raw_file_size, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.tissue)
