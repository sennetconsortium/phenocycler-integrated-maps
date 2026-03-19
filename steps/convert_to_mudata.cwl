cwlVersion: v1.1
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

requirements:
  DockerRequirement:
    dockerPull: sennet/phenocycler-maps
baseCommand: /opt/convert_to_mudata.py

inputs:
    processed_h5ad:
        type: File
        inputBinding:
            position: 0

    updated_metadata_json:
        label: "metadata about the map"
        type: File
        inputBinding:
            position: 1

    tissue:
        label: "tissue type"
        type: string
        inputBinding:
            position: 2

outputs:
    processed_h5ad:
        type: File
        outputBinding:
            glob: "*_processed.h5ad"
        doc: h5ad file with processed phenocycler datasets

    final_metadata_json:
        type: File
        outputBinding:
            glob: "*.json"
        doc: json containing data product info