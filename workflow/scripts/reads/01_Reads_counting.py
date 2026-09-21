################################################################################
# Project : "MicrobExplorer"
# Script: "Counting raw reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging
import pandas as pd

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def kaiju_analyze(PATH_IN: str, PATH_OUT: str) -> bool:
    """Parses raw Kaiju output, filters classified reads, and counts taxon occurrences."""
    try:
        df = pd.read_csv(
            PATH_IN,
            sep="\t",
            header=None,
            usecols=[0, 2],
            names=["status", "taxon_id"],
            dtype={"status": str, "taxon_id": str},
        )

        # Filtering ('C')
        df_classified = df[df["status"] == "C"]

        if not df_classified.empty:
            counts = df_classified["taxon_id"].value_counts().reset_index()

            counts.columns = ["read_id", "count"]

            counts.to_csv(PATH_OUT, sep="\t", index=False)
            return True

        return False

    except Exception as e:
        print(f"❌ Error processing Kaiju file {PATH_IN}: {e}")
        return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.raw_data
    PATH_OUT = snakemake.output.counted

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = kaiju_analyze(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[READS_COUNTING] SUCCESS | Sample: {sample_name} | "
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_COUNTING] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")