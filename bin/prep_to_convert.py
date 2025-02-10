#!/usr/bin/env python3

import os, shutil
from argparse import ArgumentParser
from pathlib import Path
import anndata


def main(processed_h5ad_file: Path):
    adata = anndata.read_h5ad(processed_h5ad_file)
    # CWL can't handle spaces in file names, so replace spaces before converting to zarr
    unds_names = []
    if "predicted_label" in adata.obs_keys():
        for name in adata.obs["predicted_label"].cat.categories:
            unds_names.append(name.replace(" ", "_"))
        adata.rename_categories("predicted_label", unds_names)

    if "cell_type_counts" in adata.obs_keys():
        cell_types = []
        for cell_type in adata.uns["cell_type_counts"]:
            cell_types.append(cell_type)
        for cell_type in cell_types:
            adata.uns["cell_type_counts"][cell_type.replace(" ", "_")] = adata.uns[
                "cell_type_counts"
            ].pop(cell_type)
    # Save where required for use with anndata-to-ui
    os.mkdir("./processed_dir")
    adata.write_h5ad("./processed_dir/secondary_analysis.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("--processed_h5ad_file", type=Path)
    args = p.parse_args()

    main(args.processed_h5ad_file)
