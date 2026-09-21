################################################################################
# Project : "MicrobExplorer"
# Script: "CPM calculation for reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def calculate_cpm(PATH_IN: str, PATH_OUT: str) -> int | bool :
    """Calculates Counts Per Million (CPM) for each taxon using Pandas."""
    try:
        df = pd.read_csv(PATH_IN, sep="\t", header=0)
        
        total_reads = df["count"].sum()
        
        if total_reads > 0:
            df["cpm"] = (df["count"] / total_reads) * 1_000_000
            
            df.to_csv(PATH_OUT, sep="\t", index=False)
            return total_reads
        return False
        
    except Exception as e:
        print(f"❌ Error processing file {PATH_IN}: {e}")
        return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.cpm

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = calculate_cpm(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[READS_CPM] SUCCESS | Sample: {sample_name} | "
            f"Total Reads for CPM: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_CPM] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"CPM calculation failed for {sample_name}")