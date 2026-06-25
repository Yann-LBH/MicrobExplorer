################################################################################
# Project : "MicrobExplorer"
# Script: "Intersection between counted contigs and kegg number extracted"
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


def intersection_kegg(PATH_IN: str, COUNTS: str, PATH_OUT: str) -> int:
    """
    Inner join between KEGG annotations and the counts for a sample.
    Returns the number of rows in the intersection.
    """
    df_kegg = pd.read_csv(PATH_IN, sep="\t")
    df_kegg["contig"] = df_kegg["contig"].astype(str).str.strip()

    df_counts = pd.read_csv(COUNTS, sep="\t")
    df_counts["contig_id"] = df_counts["contig_id"].astype(str).str.strip()
    df_counts.drop(columns=["reads_unmapped"], errors="ignore", inplace=True)

    if "Contig_ID" not in df_counts.columns:
        raise KeyError(f"Column 'contig_id' missing in {COUNTS}")

    df_out = df_counts.merge(
        df_kegg, left_on="contig_id", right_on="contig", how="inner"
    )

    if df_out.empty:
        print(f"⚠️ Aucun contig commun pour {COUNTS}")

    df_out.to_csv(PATH_OUT, sep="\t", index=False)
    return len(df_out)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    COUNTS = snakemake.input.counted
    PATH_OUT = snakemake.output.intersec

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = intersection_kegg(PATH_IN, COUNTS, PATH_OUT)
    if process:
        logging.info(
            f"[KEGG_INTERSECT] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_INTERSECT] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
