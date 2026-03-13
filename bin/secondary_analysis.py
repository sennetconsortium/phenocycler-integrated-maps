#!/usr/bin/env python3

from argparse import ArgumentParser
from matplotlib import cm
from pathlib import Path

import anndata
import json
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
import os
import pandas as pd
import scanpy as sc


def add_file_sizes(data_product_metadata, raw_size):
    data_product_metadata["Raw File Size"] = raw_size


def add_cell_counts(integrated_map_metadata, total_cell_count):
    integrated_map_metadata["Processed Total Cell Count"] = total_cell_count
    return integrated_map_metadata


def main(
    raw_h5mu_file: Path,
    integrated_map_metadata: Path,
    tissue: str = None,
):
    processed_output_file_name = (
        f"{tissue}_processed.h5mu" if tissue else "phenocycler_processed.h5mu"
    )
    # Open files and extract necessary information
    raw_mudata = md.read(raw_h5mu_file)
    with open(integrated_map_metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Integrated Map UUID"]
    adata = raw_mudata.mod[f'{uuid}_raw']

    print("Processing integrated map...")
    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.obs["n_counts"] = adata.X.sum(axis=1)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["unscaled"] = adata.X.copy()
    sc.pp.scale(adata, max_value=10)

    sc.pp.neighbors(adata, n_neighbors=50, n_pcs=50)
    sc.tl.umap(adata)

    # leiden clustering
    sc.tl.leiden(adata)
    sc.tl.rank_genes_groups(adata, "leiden")

    total_cell_count = adata.obs.shape[0]
    metadata = add_cell_counts(
            metadata, total_cell_count
        )

    with plt.rc_context():
        sc.pl.umap(adata, color="leiden", show=False)
        plt.savefig(f"{uuid}.png", bbox_inches="tight")

    # Convert to MuData and add Obj x Analyte requirements
    if 'annotation' in adata.obsm_keys():
        adata.obsm['annotation']['leiden'] = adata.obs['leiden']
    else:
        adata.obsm['annotation'] = pd.DataFrame(adata.obs['leiden'])
    adata.obsm['leiden'] = pd.DataFrame(adata.obs['leiden'])
    adata.uns['leiden'] = {
        'label': 'Leiden Clusters',
        'mechanism': 'machine',
        'protocol': "10.1186/s13059-017-1382-0",
    }

    mdata = md.MuData({f"{uuid}_processed": adata})
    mdata.uns["epic_type "] = ['analyses', 'annotations']

    print(f"Writing {processed_output_file_name}")
    mdata.write(processed_output_file_name)
    processed_file_size = os.path.getsize(processed_output_file_name)
    add_file_sizes(metadata, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5mu_file", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("integrated_map_metadata", type=Path)
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.raw_h5mu_file, args.integrated_map_metadata, args.tissue)