################################################################################
# Project : "MicrobExplorer"
# Script: "Counting raw reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def kaiju_analyze(PATH_IN, PATH_OUT):
    counter = {}

    with open(PATH_IN, "r", encoding="utf-8") as f:
        for line in f:
            # Kaiju sépare généralement les column par des tabulations
            column = line.strip().split("\t")

            # On vérifie si la line est classée ('C') et possède bien un ID taxon
            if len(column) >= 3 and column[0] == "C":
                taxon_id = column[2]  # La 3ème colonne (ex: 35786)

                counter[taxon_id] = counter.get(taxon_id, 0) + 1

    # Écriture des résultats
    if counter:
        with open(PATH_OUT, "w", encoding="utf-8") as f_out:
            # On trie par nombre d'occurrences (du plus grand au plus petit)
            for taxon, total in sorted(
                counter.items(), key=lambda x: x[1], reverse=True
            ):
                f_out.write(f"Taxon {taxon} : {total} reads\n")
        return True
    return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.raw_data
    PATH_OUT = snakemake.output.counted

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = kaiju_analyze(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[READS_COUNTING] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_COUNTING] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
