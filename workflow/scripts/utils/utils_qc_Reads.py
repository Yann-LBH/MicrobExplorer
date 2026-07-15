################################################################################
# Project : "MicrobExplorer"
# Script: "Utils for Reads QC"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os

def get_total_reads(path):
    """Count non-empty lines in a TSV or raw file, skipping the header if present."""
    count = 0
    try:
        with open(path, 'r', encoding='utf-8') as f:
            # Check if file has content and skip header if it looks like a table
            first_line = f.readline()
            if not first_line:
                return 0
                
            # If the first line contains common header keywords, treat it as header
            if not ("classified" in first_line or "unclassified" in first_line or first_line.startswith("C\t") or first_line.startswith("U\t")):
                # It's likely a standard TSV header row, so we start counting from next lines
                pass
            else:
                # It's a raw data line (e.g. Kaiju), count it
                if first_line.strip():
                    count += 1

            for line in f:
                if line.strip():
                    count += 1
    except Exception as e:
        print(f"Error reading file {path}: {e}")
        return 0
    return count


def run_full_qc(snakemake, steps_config):
    """Main execution function for Reads QC using Snakefile configuration."""
    DATA = snakemake.input.data
    qc_results = {}

    for path in DATA:
        filename = os.path.basename(path)
        normalized_path = os.path.normpath(path)

        # 1. Robust Step Detection using normalized paths
        step_key = "Unknown"
        for key, (folder, pattern) in steps_config.items():
            normalized_folder = os.path.normpath(folder)
            if normalized_folder in normalized_path:
                step_key = key
                break

        # 2. Clean sample name (Only handling TSV and Kaiju raw inputs)
        sample_name = filename
        for item in [
            "reads_",
            "counted_",
            "filtered_",
            "cpm_",
            "annotated_",
            "_reads",
            ".kaijuNR",
            ".tsv"
        ]:
            sample_name = sample_name.replace(item, "")

        # 3. Initialize sample sub-dictionary if not present
        if sample_name not in qc_results:
            qc_results[sample_name] = {}

        # 4. Force TSV/raw format parsing by setting is_csv=False
        qc_results[sample_name][step_key] = get_total_reads(path)

    return qc_results
