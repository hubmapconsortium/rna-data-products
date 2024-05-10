#!/usr/bin/env python3

import json
from argparse import ArgumentParser
from collections import defaultdict
from datetime import datetime
from os import fspath, walk, listdir
from pathlib import Path
from typing import Dict, Tuple

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import scanpy as sc
import scipy.sparse
import yaml

GENE_MAPPING_DIRECTORIES = [
    Path(__file__).parent.parent / "data",
    Path("/opt/data"),
]


annotation_fields = [
    "azimuth_label",
    "azimuth_id",
    "predicted_CLID",
    "predicted_label",
    "cl_match_type",
    "prediction_score",
]


def get_tissue_type(dataset: str) -> str:
    organ_dict = yaml.load(open("/opt/organ_types.yaml"), Loader=yaml.BaseLoader)
    organ_code = requests.get(
        f"https://entity.api.hubmapconsortium.org/dataset/{dataset}/organs/"
    )
    organ_name = organ_dict[organ_code]
    return organ_name.replace(" (Left)", "").replace(" (Right)", "")


def get_inverted_gene_dict():
    inverted_dict = defaultdict(list)
    gene_mapping = read_gene_mapping()
    for ensembl, hugo in gene_mapping.items():
        inverted_dict[hugo].append(ensembl)
    return inverted_dict


def find_files(directory, patterns):
    for dirpath_str, dirnames, filenames in walk(directory):
        dirpath = Path(dirpath_str)
        for filename in filenames:
            filepath = dirpath / filename
            for pattern in patterns:
                if filepath.match(pattern):
                    return filepath


def find_file_pairs(directory):
    filtered_patterns = ["cluster_marker_genes.h5ad", "secondary_analysis.h5ad"]
    unfiltered_patterns = ["out.h5ad", "expr.h5ad"]
    filtered_file = find_files(directory, filtered_patterns)
    unfiltered_file = find_files(directory, unfiltered_patterns)
    return filtered_file, unfiltered_file


def get_dataset_cluster_and_cell_type_if_present(barcode, filtered_adata, dataset_uuid):
    annotation_dict = {
        annotation_field: np.nan for annotation_field in annotation_fields
    }
    annotation_dict["dataset_leiden"] = np.nan
    if barcode not in filtered_adata.obs.index:
        return annotation_dict
    else:
        annotation_dict[
            "dataset_leiden"
        ] = f"{dataset_uuid}-{filtered_adata.obs.at[barcode, 'leiden']}"
        for field in annotation_fields:
            if field not in filtered_adata.obs_keys():
                annotation_dict[field] = np.nan
            else:
                annotation_dict[field] = filtered_adata.obs.at[barcode, field]
        return annotation_dict


def annotate_file(
    filtered_file: Path, unfiltered_file: Path, tissue_type: str, uuids_df: pd.DataFrame
) -> Tuple[anndata.AnnData, anndata.AnnData]:
    # Get the directory
    data_set_dir = fspath(unfiltered_file.parent.stem)
    # And the tissue type
    tissue_type = tissue_type if tissue_type else get_tissue_type(data_set_dir)
    hubmap_id = uuids_df.loc[uuids_df["uuid"] == data_set_dir, "hubmap_id"].values[0]
    filtered_adata = anndata.read_h5ad(filtered_file)
    unfiltered_adata = anndata.read_h5ad(unfiltered_file)
    unfiltered_copy = unfiltered_adata.copy()
    unfiltered_copy.obs["barcode"] = unfiltered_adata.obs.index
    unfiltered_copy.obs["dataset"] = data_set_dir
    unfiltered_copy.obs["hubmap_id"] = hubmap_id
    unfiltered_copy.obs["organ"] = tissue_type
    unfiltered_copy.obs["modality"] = "rna"
    unfiltered_copy.uns["annotation_metadata"] = (
        filtered_adata.uns["annotation_metadata"]
        if "annotation_metadata" in filtered_adata.uns.keys()
        else {"is_annotated": False}
    )
    
    # Only add annotation fields if they exist in the filtered AnnData object
    for field in annotation_fields:
        if field in filtered_adata.obs_keys():
            unfiltered_copy.obs[field] = pd.Series(
                index=unfiltered_copy.obs.index, dtype=str
            )
    if "prediction_score" in filtered_adata.obs_keys():
            unfiltered_copy.obs["prediction_score"] = pd.Series(
                index=unfiltered_copy.obs.index, dtype=np.float64
            )
            unfiltered_copy.obs["prediction_score"] = unfiltered_copy.obs["prediction_score"].astype(float)

    if any(field in unfiltered_copy.obs_keys() for field in annotation_fields):
        unfiltered_copy.obs["dataset_leiden"] = pd.Series(
            index=unfiltered_copy.obs.index, dtype=str
        )
        for barcode in unfiltered_copy.obs.index:
            dataset_clusters_and_cell_types = get_dataset_cluster_and_cell_type_if_present(
                barcode, filtered_adata, data_set_dir
            )
            for k in dataset_clusters_and_cell_types:
                unfiltered_copy.obs.at[barcode, k] = dataset_clusters_and_cell_types[k]


    cell_ids_list = [
        "-".join([data_set_dir, barcode]) for barcode in unfiltered_copy.obs["barcode"]
    ]
    unfiltered_copy.obs["cell_id"] = pd.Series(
        cell_ids_list, index=unfiltered_copy.obs.index, dtype=str
    )
    unfiltered_copy.obs.set_index("cell_id", drop=True, inplace=True)
    unfiltered_copy = map_gene_ids(unfiltered_copy)
    return unfiltered_copy



def read_gene_mapping() -> Dict[str, str]:
    """
    Try to find the Ensembl to HUGO symbol mapping, with paths suitable
    for running this script inside and outside a Docker container.
    :return:
    """
    for directory in GENE_MAPPING_DIRECTORIES:
        mapping_file = directory / "ensembl_to_symbol.json"
        if mapping_file.is_file():
            with open(mapping_file) as f:
                return json.load(f)
    message_pieces = ["Couldn't find Ensembl → HUGO mapping file. Tried:"]
    message_pieces.extend(f"\t{path}" for path in GENE_MAPPING_DIRECTORIES)
    raise ValueError("\n".join(message_pieces))


def map_gene_ids(adata):
    obsm = adata.obsm
    uns = adata.uns
    gene_mapping = read_gene_mapping()
    has_hugo_symbol = [gene in gene_mapping for gene in adata.var.index]
    # adata = adata[:, has_hugo_symbol]
    temp_df = pd.DataFrame(
        adata.X.todense(), index=adata.obs.index, columns=adata.var.index
    )
    aggregated = temp_df.groupby(level=0, axis=1).sum()
    adata = anndata.AnnData(aggregated, obs=adata.obs)
    adata.var["hugo_symbol"] = [
        gene_mapping.get(var, np.nan) for var in adata.var.index
    ]
    adata.obsm = obsm
    adata.uns = uns
    # This introduces duplicate gene names, use Pandas for aggregation
    # since anndata doesn't have that functionality
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.var_names_make_unique()
    return adata


def main(data_directory: Path, uuids_file: Path, tissue: str = None):
    raw_output_file_name = f"{tissue}_raw.h5ad" if tissue else "rna_raw.h5ad"
    processed_output_file_name = (
        f"{tissue}_processed.h5ad" if tissue else "rna_processed.h5ad"
    )
    uuids_df = pd.read_csv(uuids_file, sep="\t", dtype=str)
    directories = [data_directory / Path(uuid) for uuid in uuids_df["uuid"]]
    # Load files
    file_pairs = [find_file_pairs(directory) for directory in directories if len(listdir(directory))>1]
    print("Annotating objects")
    adatas = [
        annotate_file(file_pair[0], file_pair[1], tissue, uuids_df)
        for file_pair in file_pairs
    ]
    annotation_metadata = {
        adata.obs.dataset.iloc[0]: adata.uns["annotation_metadata"] for adata in adatas
    }
    saved_var = adatas[0].var
    print("Concatenating objects")
    adata = anndata.concat(adatas, join="outer")
    adata.uns["annotation_metadata"] = annotation_metadata
    adata.uns["creation_date_time"] = str(datetime.now())
    adata.uns["datasets"] = list(set(adata.obs.hubmap_id))
    if "predicted_label" in adata.obs_keys():
        adata.uns["cell_type_counts"] = (adata.obs["predicted_label"].value_counts().to_dict())
    adata.var = saved_var
    print(f"Writing {raw_output_file_name}")
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
    p.add_argument("data_directory", type=Path)
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.data_directory, args.uuids_file, args.tissue)