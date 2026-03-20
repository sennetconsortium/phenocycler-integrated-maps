#!/usr/bin/env python3

from argparse import ArgumentParser
from matplotlib import cm
from pathlib import Path

import anndata as ad
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import rapids_singlecell as rsc
import scanpy as sc


def add_cell_counts(integrated_map_metadata, total_cell_count):
    integrated_map_metadata["Processed Total Cell Count"] = total_cell_count
    return integrated_map_metadata


def main(
    raw_h5ad_file: Path,
    integrated_map_metadata: Path,
    tissue: str = None,
):
    processed_output_file_name = (
        f"{tissue}_processed" if tissue else "phenocycler_processed"
    )
    # Open files and extract necessary information
    adata = ad.read_h5ad(raw_h5ad_file)
    with open(integrated_map_metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Integrated Map UUID"]

    # Move .X to the GPU for rsc
    rsc.get.anndata_to_GPU(adata)

    print("Processing integrated map...")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    adata.obs["n_counts"] = adata.X.sum(axis=1).get()

    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)
    adata.layers["unscaled"] = adata.X.copy()
    rsc.pp.scale(adata, max_value=10)

    rsc.pp.neighbors(adata, n_neighbors=50)
    rsc.tl.umap(adata)

    # Move .X back to the CPU to use regular sc for leiden
    rsc.get.anndata_to_CPU(adata)
    # leiden clustering, maybe try the use dask param with rapids?
    sc.tl.leiden(adata, flavor='igraph')

    total_cell_count = adata.obs.shape[0]
    metadata = add_cell_counts(
            metadata, total_cell_count
        )
    with open(f"{uuid}.json", "w") as outfile:
        json.dump(metadata, outfile)

    # Move .X back to the CPU to plot
    # rsc.get.anndata_to_CPU(adata)

    # Plot
    with plt.rc_context():
        sc.pl.umap(adata, color='leiden', show=False)
        plt.savefig(f"{uuid}.png", bbox_inches="tight")

    # # Convert to MuData and add Obj x Analyte requirements
    # if 'annotation' in adata.obsm_keys():
    #     adata.obsm['annotation']['leiden'] = adata.obs['leiden']
    # else:
    #     adata.obsm['annotation'] = pd.DataFrame(adata.obs['leiden'])
    # adata.obsm['leiden'] = pd.DataFrame(adata.obs['leiden'])
    # adata.uns['leiden'] = {
    #     'label': 'Leiden Clusters',
    #     'mechanism': 'machine',
    #     'protocol': "10.1186/s13059-017-1382-0",
    # }

    print(f"Writing {processed_output_file_name}")
    adata.write(f"{processed_output_file_name}.h5ad")


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5mu_file", type=Path)
    p.add_argument("integrated_map_metadata", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.raw_h5mu_file, args.integrated_map_metadata, args.tissue)