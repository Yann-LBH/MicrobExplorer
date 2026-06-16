################################################################################
# Project : "MicrobExplorer"
# Script: "Calcul RPKM for contigs"
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


def calculate_rpkm(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Calcul the contigs RPKM
    Formula : (Reads * 10^9) / (Length * Total_Mapped_Reads)
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    required_cols = {"Reads_Mapped", "Length"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes dans {PATH_IN} : {missing}")

    total_mapped = df["Reads_Mapped"].sum()

    if total_mapped == 0:
        print(f"⚠️  Total mapped reads = 0 pour {PATH_IN}, RPKM mis à NaN.")
        df["RPKM"] = float("nan")
        df.to_csv(PATH_OUT, sep="\t", index=False)
        return 0

    df["RPKM"] = ((df["Reads_Mapped"] * 1e9) / (df["Length"] * total_mapped)).round(4)
    df.to_csv(PATH_OUT, sep="\t", index=False)
    return len(df)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.rpkm

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = calculate_rpkm(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[CONTIGS_RPKM] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[CONTIGS_RPKM] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
