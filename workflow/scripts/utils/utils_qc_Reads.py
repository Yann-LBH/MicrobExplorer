################################################################################
# Project : "MicrobExplorer"
# Script: "Utils for Reads QC"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################


def run_full_qc(snakemake, steps_config):
    """Main execution function for Reads QC using Snakefile configuration."""
    DATA = snakemake.input.data
    qc_results = {}

    # steps_config looks like: {"brut": ("Data", "reads_{sample}.kaijuNR"), ...}
    for path in DATA:
        filename = os.path.basename(path)

        # Detect the step by matching the folder name or pattern from your Snakefile dictionary
        step_key = "Unknown"
        for key, (folder, pattern) in steps_config.items():
            if folder in path:
                step_key = key
                break

        # Clean sample name
        sample_name = filename
        # Basic cleanup: remove extensions and prefixes to isolate the sample name
        for item in [
            "counted_",
            "filtered_",
            "annotated_",
            "reads_",
            "_reads",
            ".tsv",
            ".csv",
            ".kaijuNR",
        ]:
            sample_name = sample_name.replace(item, "")

        # Call your get_total_reads logic...
        is_csv = step_key == "final"  # True if it's the final annotated table

        if sample_name not in qc_results:
            qc_results[sample_name] = {}

        qc_results[sample_name][step_key] = get_total_reads(path, is_csv=is_csv)

    return qc_results
