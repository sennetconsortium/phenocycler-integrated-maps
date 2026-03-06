#!/usr/bin/env python3
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd
import re
import requests
import subprocess


def download_from_assets(uuid, ancestor_uuid, file_patterns, out_dir, bearer_token = None):
    files_directory = out_dir / uuid
    files_directory.mkdir(parents=True, exist_ok=True)

    if bearer_token:
        headers = {"Authorization": f"Bearer {bearer_token}"}
    else:
        headers = None

    for pattern in file_patterns:
        print(pattern)
        out_file = files_directory / pattern
        url = f"https://assets.api.sennetconsortium.org/{uuid}/sprm_outputs/{pattern}"
        # Check if file exists and download
        try:
            head = requests.head(url, allow_redirects=True, headers=headers)
            if head.status_code != 200:
                print(f"{pattern} for {uuid} does not exist on server (status {head.status_code}).")
                continue
        except requests.RequestException as e:
            print(f"Error checking {uuid}: {e}")
            continue
        try:
            print(f"Downloading {pattern} for {uuid}...")
            with requests.get(url, stream=True, headers=headers) as r:
                r.raise_for_status()
                with open(out_file, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        f.write(chunk)
            print(f"Downloaded {pattern} for {uuid} successfully.")
        except requests.RequestException as e:
            print(f"Failed to download {uuid}: {e}")

    # TODO: Antibodies tsv time, need to figure out regex instead of hardcode
    antibodies_url = f"https://assets.api.sennetconsortium.org/{ancestor_uuid}/extras/SenNet_antibodies_n40_060325.tsv"
    out = Path(files_directory / "SenNet_antibodies_n40_060325.tsv")
    try:
        head = requests.head(antibodies_url, allow_redirects=True, headers=headers)
        if head.status_code != 200:
            print(f"Couldn't get antibodies for processed dataset {uuid}, might be filename issue.")
            print(f"Tried path {antibodies_url}")
    except requests.RequestException as e:
        print(f"Error getting antibodies for processed dataset {uuid}: {e}")
        print(f"Tried path {antibodies_url}")
    try:
        print(f"Found antibodies for processed dataset {uuid}, downloading {antibodies_url}...")
        with requests.get(antibodies_url, stream=True, headers=headers) as r:
            r.raise_for_status()
            with open(out, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
                print(f"Downloaded antibodies tsv successfully.")
    except requests.RequestException as e:
        print(f"Failed to download antibodies for processed dataset {uuid}: {e}")


# def find_antibodies_tsv(uuid, pattern, token = None):
#     if token:
#         headers = {"Authorization": f"Bearer {token}"}
#     else:
#         headers = None
#     url = f"https://assets.api.sennetconsortium.org/{uuid}/{pattern}"
#     response = requests.get(url, allow_redirects=True, stream=True, headers=headers)


def main(uuids_file: Path, tissue: str, token = None):
    df = pd.read_csv(uuids_file, sep = "\t")
    processed_uuids = df["uuid"]
    ancestor_uuids = df["ancestors"]
    files_base_directory = Path(f"{tissue}_files")
    files_base_directory.mkdir(exist_ok=True)

    # SPRM outputs first, then antibodies tsv
    hdf5_pattern = Path("out.hdf5")
    cell_count_pattern = Path("aligned_tissue_0_expr.ome.tiff-cell_channel_total.csv")
    adjacency_matrix_pattern = Path("aligned_tissue_0_expr.ome.tiff_AdjacencyMatrix.mtx")
    adjacency_matrix_labels_pattern = Path("aligned_tissue_0_expr.ome.tiff_AdjacencyMatrixRowColLabels.txt")
    cell_centers_pattern = Path("aligned_tissue_0_expr.ome.tiff-cell_centers.csv")
    # antibodies_re = re.compile(r"^extras/.*antibodies.*\.tsv$")
    file_patterns = [hdf5_pattern, cell_count_pattern, adjacency_matrix_pattern, adjacency_matrix_labels_pattern, cell_centers_pattern]
    for uuid, ancestor in zip(processed_uuids, ancestor_uuids):
        # TODO: need to figure out how to find the antibodies tsv with just a regex and not hardcode
        # antibodies_pattern = find_antibodies_tsv(ancestor, antibodies_re, token)
        # file_patterns.append("extras/SenNet_antibodies_n40_060325.tsv")
        download_from_assets(uuid, ancestor, file_patterns, files_base_directory, token)



if __name__ == "__main__":
    p = ArgumentParser()
    p.add_argument("uuids_file", type=Path)
    p.add_argument("tissue", type=str)
    p.add_argument("bearer_token", type=str, nargs="?")

    args = p.parse_args()

    main(args.uuids_file, args.tissue, args.bearer_token)