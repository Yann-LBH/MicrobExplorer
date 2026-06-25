################################################################################
# Project : "MicrobExplorer"
# Script: "Data prepared for DEseq2"
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

METADATA_COLS = {"contig", "length", "contig_id", "gene_length", "id", "product"}


def sum_per_kegg(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Aggregates the counts by KEGG ID.
    Returns the number of unique KOs.
    """
    header = pd.read_csv(PATH_IN, sep="\t", nrows=0).columns
    cols_keep = [c for c in header if c == "kegg" or c not in METADATA_COLS]

    df = pd.read_csv(PATH_IN, sep="\t", usecols=cols_keep)

    if "kegg" not in df.columns:
        raise KeyError(f"Column 'kegg' missing in {PATH_IN}")

    df_out = df.groupby("kegg").sum(numeric_only=True).reset_index()
    df_out.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_out)


# --- Exécution ---
if __name__ == "__main__":
    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.deseq2

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = sum_per_kegg(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[KEGG_DESEQ2] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_DESEQ2] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
