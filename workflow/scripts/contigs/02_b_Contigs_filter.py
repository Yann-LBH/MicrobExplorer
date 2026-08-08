################################################################################
# Project : "MicrobExplorer"
# Script: "Filter the contigs abundance"
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
    PATH_IN = str(snakemake.input.data)
    GLOBAL_ABUNDANCE = str(snakemake.input.global_abundance)
    PATH_OUT = str(snakemake.output.filtered)
    THRESHOLD = int(snakemake.params.abundance_threshold)

    # 1. Load the pre-calculated global reference counts
    ref_df = pd.read_csv(GLOBAL_ABUNDANCE, sep="\t", index_col="contig_id")

    # 2. Load and filter the single current sample file
    df = pd.read_csv(PATH_IN, sep="\t")

    # Map the target contigs to the pre-computed global values
    df_filtered = df[
        df["contig_id"].map(ref_df["total_abundance"]).fillna(0) >= THRESHOLD
    ]

    # 3. Save the individual filtered file
    try:
        df_filtered.to_csv(PATH_OUT, sep="\t", index=False)
        logging.info(
            f"[CONTIGS_FILTER] SUCCESS | Output: {PATH_OUT}"
        )
    except Exception as e:
        logging.error(
            f"[CONTIGS_FILTER] FAILED | Output: {PATH_OUT} | Error: {e}"
        )
        raise
