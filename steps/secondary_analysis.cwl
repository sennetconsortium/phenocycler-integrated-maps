cwlVersion: v1.0
class: CommandLineTool
label: Annotates each h5ad file with dataset and tissue type, then concatenates

hints:
  DockerRequirement:
    dockerPull: sennet/phenocycler-maps
baseCommand: /opt/concatenate.py

inputs:
    raw_h5mu:
        label: "Where the h5ad files are"
        type: File
        inputBinding:
            position: 0

    metadata_json:
        label: "metadata about the map"
        type: File
        inputBinding:
            position: 0


outputs:
    processed_h5mu:
        type: File
        outputBinding:
            glob: "*_processed.h5mu"
        doc: h5mu file with concatenated codex datasets

    final_metadata_json:
        type: File
        outputBinding:
            glob: "*.json"
        doc: json containing data product info