################################################################################
# Project : "MicrobExplorer"
# Script: "Counting contigs"
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


# --- Filtering ---
def filter_contigs(PATH_IN: str, PATH_OUT: str, THRESHOLD: int) -> int:
    """
    Reads the input, filters based on `Length` >= `threshold`,
    and writes the output to a TSV file.
    Returns the number of retained contigs.
    """
    df = pd.read_csv(
        PATH_IN, 
        sep="\t", 
        header=None, 
        names=["Contig_id", "Length", "Reads_Mapped", "Reads_Unmapped"],
        dtype={"Length": int}
    )

    df_filtered = df[df["Length"] >= THRESHOLD]

    df_filtered.to_csv(PATH_OUT, sep="\t", index=False, header=True)
    return len(df_filtered)


# --- Execution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.raw_data
    PATH_OUT = snakemake.output.counted
    THRESHOLD = int(snakemake.params.length_threshold)

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = filter_contigs(PATH_IN, PATH_OUT, THRESHOLD)
    if process:
        logging.info(
            f"[CONTIGS_COUNTING] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[CONTIGS_COUNTING] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
