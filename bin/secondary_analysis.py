#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path

import anndata
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def add_cell_counts(data_product_metadata, cell_counts, total_cell_count):
    with open(data_product_metadata, 'r') as json_file:
        metadata = json.load(json_file)
    uuid = metadata["Data Product UUID"]
    metadata["Processed Cell Type Counts"] = cell_counts
    metadata["Processed Total Cell Count"] = total_cell_count
    with open(f"{uuid}.json", 'w') as outfile:
        json.dump(metadata, outfile)


def main(raw_h5ad_file: Path, data_product_metadata: Path, tissue: str = None):
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )

    adata = anndata.read_h5ad(raw_h5ad_file)
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
        add_cell_counts(data_product_metadata, cell_type_counts, total_cell_count)
    else:
        cell_type_counts = {}
        add_cell_counts(data_product_metadata, cell_type_counts, total_cell_count)

    with plt.rc_context():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig(f"{uuid}.png" if tissue else "rna.png")

    print(f"Writing {processed_output_file_name}")
    adata.write(processed_output_file_name)


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