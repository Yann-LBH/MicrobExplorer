################################################################################
# Project : "MicrobExplorer"
# Script: "Calcule the contigs abundance"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import logging
import pandas as pd

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

if __name__ == "__main__":
    PATH_IN = snakemake.input.all_data
    PATH_OUT = str(snakemake.output.global_abundance)

    global_counts = pd.Series(dtype=int)

    # Sum abundances across all files
    for f in PATH_IN:
        chunk = pd.read_csv(f, sep="\t", usecols=["contig_id", "read_mapped"])

        file_counts = chunk.groupby("contig_id")["read_mapped"].sum()
        global_counts = global_counts.add(file_counts, fill_value=0)

    # Save the reference global table
    try:
        global_counts.to_csv(
            PATH_OUT, sep="\t", header=["total_abundance"], index_label="contig_id"
        )
        logging.info(
            f"[CONTIGS_GLOBAL_ABUNDANCE] SUCCESS | Output: {PATH_OUT}"
        )
    except Exception as e:
        logging.error(
            f"[CONTIGS_GLOBAL_ABUNDANCE] FAILED | Output: {PATH_OUT} | Error: {e}"
        )
        raise
