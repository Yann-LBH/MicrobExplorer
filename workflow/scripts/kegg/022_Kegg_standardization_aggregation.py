################################################################################
# Project : "MicrobExplorer"
# Script: "Standardization of aggregation by unique KO per files"
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


def aggregate_by_ko(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Aggregate by KEGG code by summing the ‘standardization’ column.
    Returns the number of unique KOs.
    """
    df = pd.read_csv(PATH_IN, sep="\t", usecols=["kegg", "standardization"])

    missing = {"kegg", "standardization"} - set(df.columns)
    if missing:
        raise KeyError(f"Column missing : {missing}")

    df_ko = df.groupby("kegg", as_index=False)["standardization"].sum()
    df_ko.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_ko)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.agreg

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = aggregate_by_ko(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[KEGG_AGGREGATION] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_AGGREGATION] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
