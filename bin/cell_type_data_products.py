#!/usr/bin/env python3

from argparse import ArgumentParser
from pathlib import Path
from typing import List

import anndata
import scanpy as sc

def main(raw_products: List[Path]):

    # Load files
    print("Annotating objects")
    adatas = [anndata.read(product) for product in raw_products]
    adata = anndata.concat(adatas)

    #Use mapping magic to add different levels of granularity of cell type assignments
    #And repeat this process at each different level

    print(adata.obs.columns)

    for cell_type in adata.obs["predicted_label"].dropna().unique():
        print(cell_type)
        cell_type_adata = adata[adata.obs["predicted_label"] == cell_type]
        if len(cell_type_adata.obs.index) < 200:
            pass

        cell_type_adata.write(f"{cell_type}_raw.h5ad")
        print("Processing data product")
        cell_type_adata.var_names_make_unique()
        cell_type_adata.obs_names_make_unique()

        sc.pp.filter_cells(cell_type_adata, min_genes=200)
        sc.pp.filter_genes(cell_type_adata, min_cells=3)

        cell_type_adata.obs["n_counts"] = cell_type_adata.X.sum(axis=1)

        sc.pp.normalize_total(cell_type_adata, target_sum=1e4)
        sc.pp.log1p(cell_type_adata)
        cell_type_adata.layers["unscaled"] = cell_type_adata.X.copy()
        #    sc.pp.combat(adata, "dataset")
        sc.pp.scale(cell_type_adata, max_value=10)

        sc.pp.pca(cell_type_adata, n_comps=50)
        sc.pp.neighbors(cell_type_adata, n_neighbors=50, n_pcs=50)

        sc.tl.umap(cell_type_adata)

        # leiden clustering
        sc.tl.leiden(cell_type_adata)


        counts_dict = cell_type_adata.obs.organ.value_counts().to_dict()
        keep_organs = [
            organ for organ in counts_dict if counts_dict[organ] > 1
        ]
        adata_filter = cell_type_adata[cell_type_adata.obs.organ.isin(keep_organs)]
        # Filter out cell types with only one cell for this analysis
        sc.tl.rank_genes_groups(adata_filter, "organ", key_added="rank_genes_groups_organ")
        cell_type_adata.uns = adata_filter.uns
        sc.tl.rank_genes_groups(adata_filter, "leiden", key_added="rank_genes_groups_leiden")

        sc.pl.umap(
            cell_type_adata, color="leiden", show=False, save=f"{cell_type}_by_leiden.png"
        )

        sc.pl.umap(
            cell_type_adata, color="organ", show=False, save=f"{cell_type}_by_organ.png"
        )

        processed_output_file_name = f"{cell_type}_processed.h5ad"
        print(f"Writing {processed_output_file_name}")
        cell_type_adata.write(processed_output_file_name)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("organ_raw_products", type=Path, nargs="+")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole

        manhole.install(activate_on="USR1")

    main(args.organ_raw_products)