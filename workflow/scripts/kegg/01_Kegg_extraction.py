################################################################################
# Project : "MicrobExplorer"
# Script: "Extraction of lines that have a KEGG number"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import re
import os
import logging

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# Pré-compilation regex — justifié car appelée sur chaque ligne
KEGG_PATTERN = re.compile(r"K\d{5}")
HEADER = ["contig", "kegg", "gene_length"]


def process_gff_kegg(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Parses a GFF3 file and extracts KEGG annotations in real time.
    Returns the number of annotations extracted.
    """
    count = 0

    with open(PATH_IN, "r") as f_in, open(PATH_OUT, "w", newline="") as f_out:

        f_out.write("\t".join(HEADER) + "\n")

        for line in f_in:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue

            try:
                length = int(parts[4]) - int(parts[3]) + 1
            except ValueError:
                continue

            for kegg in set(KEGG_PATTERN.findall(parts[8])):
                f_out.write(f"{parts[0]}\t{kegg}\t{length}\n")
                count += 1

    return count


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.gff
    PATH_OUT = snakemake.output.tsv

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = process_gff_kegg(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[KEGG_EXTRACT] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_EXTRACT] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
