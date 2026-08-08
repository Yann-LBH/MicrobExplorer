################################################################################
# Project : "MicrobExplorer"
# Script: "Filter contigs by RPKM treshold"
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


def get_min_rpkm_across_samples(files: list[str]) -> pd.Series:
    """
    Returns the minimum RPKM per contig across all samples.
    A contig that is absent from a sample is assigned a value of 0.0 for that sample.
    """
    frames = [
        pd.read_csv(f, sep="\t", usecols=["contig_id", "rpkm"]).set_index("contig_id")[
            "rpkm"
        ]
        for f in files
    ]
    # Create a table with sample in col and contig_id + rpkm in row. 
    # Take the minimal rpkm value of all files and fill the row with.
    return pd.concat(frames, axis=1).fillna(0.0).min(axis=1)


def filter_by_min_rpkm(
    current: str, PATH_OUT: str, global_min: pd.Series, RPKM_THRESHOLD: float
) -> int:
    """
    Filters the current TSV file: keeps contigs with a minimum global RPKM value >= to the threshold.
    Returns the number of lines retained.
    """
    df = pd.read_csv(current, sep="\t")
    # Mask True / False based on global value
    mask = df["contig_id"].map(global_min).fillna(0.0) >= RPKM_THRESHOLD
    # Filter original df with mask
    df[mask].to_csv(PATH_OUT, sep="\t", index=False)
    return mask.sum()


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input
    PATH_OUT = str(snakemake.output.rpkm_filtered)
    RPKM_THRESHOLD = float(snakemake.params.rpkm_threshold)

    current = snakemake.input[0]

    global_min = get_min_rpkm_across_samples(PATH_IN)

    # Report
    sample_name = snakemake.wildcards.sample
    process = filter_by_min_rpkm(current, PATH_OUT, global_min, RPKM_THRESHOLD)
    if process:
        logging.info(
            f"[CONTIGS_RPKM_FILTER] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[CONTIGS_RPKM_FILTER] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
