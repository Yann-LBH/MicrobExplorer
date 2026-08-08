################################################################################
# Project : "MicrobExplorer"
# Script: "Filter for reads"
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


# All script comments are provided in English as requested.
def filtering_kaiju(PATH_IN: str, PATH_OUT: str, COUNT_THRESHOLD: int) -> bool:
    """Filters Kaiju output lines based on a minimum read count threshold."""
    try:
        # 1. On charge le fichier en conservant le header d'origine intact
        df = pd.read_csv(PATH_IN, sep="\t", header=0)
        
        # On cible la deuxième colonne (le count) de manière dynamique via son index
        count_col = df.columns[1]
        
        # 2. Filtrage vectorisé ultra-rapide
        df_filtered = df[df[count_col] >= COUNT_THRESHOLD]
        
        # 3. Sauvegarde propre
        df_filtered.to_csv(PATH_OUT, sep="\t", index=False)
        return True

    except Exception as e:
        print(f"❌ Error processing file {PATH_IN}: {e}")
        return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = str(snakemake.input.data)
    PATH_OUT = str(snakemake.output.filtered)
    COUNT_THRESHOLD = int(snakemake.params.count_threshold)

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = filtering_kaiju(PATH_IN, PATH_OUT, COUNT_THRESHOLD)
    if process:
        logging.info(
            f"[READS_FILTER] SUCCESS | Sample: {sample_name} | "
            f"Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_FILTER] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")