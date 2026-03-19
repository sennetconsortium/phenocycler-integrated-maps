#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path
import json
import os

import anndata as ad
import mudata as md

def add_file_sizes(data_product_metadata, raw_size):
    data_product_metadata["Processed File Size"] = raw_size


def main(processed_h5ad, metadata, tissue=None):
    output_file_name = f"{tissue}_processed" if tissue else "phenocycler_processed"
    adata = ad.read_h5ad(processed_h5ad)
    with open(metadata, "r") as infile:
        metadata = json.load(infile)
    uuid = metadata["Integrated Map UUID"]
    mdata = md.MuData({f"{uuid}_processed": adata})
    mdata.uns["epic_type "] = ['analyses', 'annotations']
    print(f"Writing {output_file_name}")
    mdata.write(f"{output_file_name}.h5mu")
    processed_file_size = os.path.getsize(f"{output_file_name}.h5mu")
    add_file_sizes(metadata, processed_file_size)


if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("raw_h5ad_file", type=Path)
    p.add_argument("integrated_map_metadata", type=Path)
    p.add_argument("tissue", type=str, nargs="?")
    p.add_argument("--enable_manhole", action="store_true")

    args = p.parse_args()

    if args.enable_manhole:
        import manhole
        manhole.install(activate_on="USR1")

    main(args.raw_h5ad_file, args.integrated_map_metadata, args.tissue)
