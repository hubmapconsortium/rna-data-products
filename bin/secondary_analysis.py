#!/usr/bin/env python3

from argparse import ArgumentParser
from os import fspath, walk, listdir
from pathlib import Path
from typing import Dict, Tuple

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


def add_azimuth_annotations(h5ad, csv):
    h5ad.obs = pd.merge(h5ad.obs, csv, "outer")
    return h5ad


def main(raw_h5ad_file, annotated_csv: None, tissue: None, metadata_json):
    raw_output_file_name = f"{tissue}_raw.h5ad" if tissue else "rna_raw.h5ad"
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )
    adata = anndata.read(raw_h5ad_file)
    if annotated_csv:
        annotations = pd.read_csv(annotated_csv)
        adata = add_azimuth_annotations(adata, annotations)
    print("Writing raw file with Azimuth annotations")
    adata.write(raw_output_file_name)    
    print("Processing data product")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

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

    # TODO: deal with predicted label

    if "predicted_label" in adata.obs_keys():

        non_na_values = adata.obs.predicted_label.dropna()
        counts_dict = non_na_values.value_counts().to_dict()
        keep_cell_types = [
            cell_type for cell_type in counts_dict if counts_dict[cell_type] > 1
        ]
        adata_filter = adata[adata.obs.predicted_label.isin(keep_cell_types)]
        # Filter out cell types with only one cell for this analysis
        sc.tl.rank_genes_groups(adata_filter, "predicted_label")
        adata.uns = adata_filter.uns
    
    if "predicted_label" in adata.obs_keys():
        adata.uns["cell_type_counts"] = (adata.obs["predicted_label"].value_counts().to_dict())

    sc.pl.umap(
        adata, color="leiden", show=False, save=f"{tissue}.png" if tissue else "rna.png"
    )
    print(f"Writing {processed_output_file_name}")
    adata.write(processed_output_file_name)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("annotated_csv", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("metadata_json", type=Path)

    args = p.parse_args()

    main(args.raw_h5ad_file, args.annotated_csv, args.tissue, args.metadata_json)