################################################################################
# Project : "MicrobExplorer"
# Script: "Union of contigs filtered by length (Script 02) and contigs filtered
#          by RPKM value ( Script 04)"
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


def get_union_and_extract(
    ABUNDANCE: str, RPKM_FILTERED: str, RAW_DATA: str, PATH_OUT: str
) -> int:
    """
    Union des Contig_ID des deux fichiers filtrés,
    puis extraction des lignes correspondantes depuis le fichier source.
    Retourne le nombre de contigs extraits.
    """
    ids_abund = pd.read_csv(ABUNDANCE, sep="\t", usecols=["Contig_ID"])["Contig_ID"]
    ids_rpkm = pd.read_csv(RPKM_FILTERED, sep="\t", usecols=["Contig_ID"])["Contig_ID"]

    target_ids = set(ids_abund).union(ids_rpkm)

    df_source = pd.read_csv(RAW_DATA, sep="\t")
    df_out = df_source[df_source["Contig_ID"].isin(target_ids)]
    df_out.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_out)


# --- Exécution ---
if __name__ == "__main__":

    ABUNDANCE = snakemake.input.abundance
    RPKM_FILTERED = snakemake.input.rpkm_filtered
    RAW_DATA = snakemake.input.raw_data
    PATH_OUT = snakemake.output.union

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = get_union_and_extract(ABUNDANCE, RPKM_FILTERED, RAW_DATA, PATH_OUT)
    if process:
        logging.info(
            f"[CONTIGS_UNION_FILTER] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[CONTIGS_UNION_FILTER] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
