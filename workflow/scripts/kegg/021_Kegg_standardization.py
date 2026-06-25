################################################################################
# Project : "MicrobExplorer"
# Script: "Standardization by weighted abundance"
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

REQUIRED_COLS = {"reads_mapped", "gene_length", "length"}


def calculate_score_kegg(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Calculates the weighted abundance: (reads_mapped * gene_length) / length
    Returns the number of processed lines.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise KeyError(f"Column missing in {PATH_IN} : {missing}")

    df["standardization"] = (df["reads_mapped"] * df["gene_length"]) / df["length"]
    df.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.stand

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = calculate_score_kegg(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[KEGG_STANDARDIZATION] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_STANDARDIZATION] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
