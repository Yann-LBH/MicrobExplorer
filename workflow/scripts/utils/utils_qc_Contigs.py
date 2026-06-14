################################################################################
# Project : "MicrobExplorer"
# Script: "Utils for Contigs QC"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################


def run_full_qc(snakemake, steps_config):
    """Main execution function for Contigs QC using Snakefile configuration."""
    DATA = snakemake.input.data
    qc_results = {}

    for path in DATA:
        filename = os.path.basename(path)

        # ✅ Dynamic detection based on your Snakefile QC_STEPS!
        step_key = "Unknown"
        for key, (folder, pattern) in steps_config.items():
            if folder in path:
                step_key = key
                break

        # Clean sample name
        sample_name = filename
        for item in [
            "count-contigs-coassembly-",
            "Raw_counting_",
            "Filtered_Contigs_",
            "RPKM_Contigs_",
            "Trimmed_",
            "Final_Selected_",
            "annotated_",
            "union_",
            "rpkm_filtered_",
            "rpkm_",
            "counted_",
            "filtered_",
            "contigs_",
            ".tsv",
            ".csv",
            ".txt",
            ".kaijuNR",
        ]:
            sample_name = sample_name.replace(item, "")

        is_table = filename.endswith(".csv") or filename.endswith(".tsv")
        ids = get_contig_ids(path, is_csv=is_table)

        if sample_name not in qc_results:
            qc_results[sample_name] = {}

        qc_results[sample_name][step_key] = len(ids)

    return qc_results
