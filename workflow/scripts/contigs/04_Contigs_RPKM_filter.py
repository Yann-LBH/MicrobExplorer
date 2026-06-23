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
    Retourne le RPKM minimum par contig sur tous les échantillons.
    Un contig absent d'un échantillon vaut 0.0 pour cet échantillon.
    """
    frames = [
        pd.read_csv(f, sep="\t", usecols=["Contig_ID", "RPKM"]).set_index("Contig_ID")[
            "RPKM"
        ]
        for f in files
    ]
    # concat en colonnes puis min par ligne — gère les NaN (contig absent) comme 0.0
    return pd.concat(frames, axis=1).fillna(0.0).min(axis=1)


def filter_by_min_rpkm(
    current: str, PATH_OUT: str, global_min: pd.Series, RPKM_THRESHOLD: float
) -> int:
    """
    Filtre le TSV courant : garde les contigs dont le RPKM min global >= threshold.
    Retourne le nombre de lignes conservées.
    """
    df = pd.read_csv(current, sep="\t")
    mask = df["Contig_ID"].map(global_min).fillna(0.0) >= RPKM_THRESHOLD
    df[mask].to_csv(PATH_OUT, sep="\t", index=False)
    return mask.sum()


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input
    PATH_OUT = snakemake.output.rpkm_filtered
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
