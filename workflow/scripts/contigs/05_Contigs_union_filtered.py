################################################################################
# Project : "MicrobExplorer"
# Script: "Union of contigs filtered by length (Script 02) and contigs filtered
#          by RPKM value ( Script 04)"
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

def get_union_and_extract(
    ABUNDANCE: str, RPKM_FILTERED: str, DATA_SOURCE: str, PATH_OUT: str
) -> int:
    """
    Union des Contig_ID des deux fichiers filtrés,
    puis extraction des lignes correspondantes depuis le fichier source.
    Retourne le nombre de contigs extraits.
    """
    
    ids_abund = pd.read_csv(ABUNDANCE, sep="\t", usecols=["contig_id"])["contig_id"]
    ids_rpkm = pd.read_csv(RPKM_FILTERED, sep="\t", usecols=["contig_id"])["contig_id"]

    target_ids = set(ids_abund).union(ids_rpkm)

    df_source = pd.read_csv(DATA_SOURCE, sep="\t")

    df_out = df_source[df_source["contig_id"].isin(target_ids)]
    df_out.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_out)


# --- Exécution ---
if __name__ == "__main__":

    ABUNDANCE = snakemake.input.abundance
    RPKM_FILTERED = snakemake.input.rpkm_filtered
    DATA_SOURCE = snakemake.input.data_source
    PATH_OUT = snakemake.output.union

    # Report
    sample_name = snakemake.wildcards.sample
    process = get_union_and_extract(ABUNDANCE, RPKM_FILTERED, DATA_SOURCE, PATH_OUT)
    if process:
        logging.info(
            f"[CONTIGS_UNION_FILTER] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[CONTIGS_UNION_FILTER] FAILED  | Sample: {sample_name}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
