################################################################################
# Project : "MicrobExplorer"
# Script: "Utils for Contigs QC"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os

def get_contig_ids(path):
    """Read a TSV file and return a list of unique contig IDs, skipping the header row."""
    ids = []
    try:
        with open(path, 'r', encoding='utf-8') as f:
            # Skip header row (e.g., 'contig_id\tcount' or similar column names)
            next(f, None)
            for line in f:
                if line.strip():
                    # Split by tab since it is exclusively a TSV file
                    parts = line.split('\t')
                    if parts[0].strip():
                        ids.append(parts[0].strip())
    except Exception as e:
        print(f"Error reading contig file {path}: {e}")
        return []
    return ids

def run_full_qc(snakemake, steps_config):
    """Main execution function for Contigs QC using Snakefile configuration."""
    DATA = snakemake.input.data
    qc_results = {}

    for path in DATA:
        filename = os.path.basename(path)
        normalized_path = os.path.normpath(path)

        # 1. Robust Step Detection
        step_key = "Unknown"
        for key, (folder, pattern) in steps_config.items():
            normalized_folder = os.path.normpath(folder)
            if normalized_folder in normalized_path:
                step_key = key
                break

        # 2. Clean sample name (Only handling TSV)
        sample_name = filename
        for item in [
            "count-contigs-coassembly-",
            "counted_",
            "filtered_",
            "rpkm_",
            "rpkm_filtered_",
            "union_",
            "annotated_",
            "_contigs",
            ".tsv"
        ]:
            sample_name = sample_name.replace(item, "")

        # 3. Force TSV format parsing
        ids = get_contig_ids(path)

        # 4. Structure results dictionary safely
        if sample_name not in qc_results:
            qc_results[sample_name] = {}

        qc_results[sample_name][step_key] = len(ids)

    return qc_results