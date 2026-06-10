################################################################################
# Project : "MicrobExplorer"
# Script: "Extraction of lines that have a KEGG number"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import re

# Pré-compilation regex — justifié car appelée sur chaque ligne
KEGG_PATTERN = re.compile(r'K\d{5}')
HEADER = ["contig", "kegg", "gene_length"]

def process_gff_kegg(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Parse un fichier GFF3 et extrait les annotations KEGG en streaming.
    Retourne le nombre d'annotations extraites.
    """
    count = 0

    with open(PATH_IN, "r") as f_in, \
         open(PATH_OUT, "w", newline="") as f_out:

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

    n = process_gff_kegg(PATH_IN, PATH_OUT)

    print(f"✓ KEGG : Extraction step passed successfully -> {n} annotations extracted -> {PATH_OUT}")