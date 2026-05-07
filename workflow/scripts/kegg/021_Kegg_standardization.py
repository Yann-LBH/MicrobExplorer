################################################################################
# Project : "MicrobExplorer"
# Script: "Standardization by weighted abundance"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################


import pandas as pd

REQUIRED_COLS = {"Reads_Mapped", "gene_length", "Length"}

def calculer_score_kegg(path_in: str, path_out: str) -> int:
    """
    Calcule l'abondance pondérée : (Reads_Mapped * gene_length) / Length
    Retourne le nombre de lignes traitées.
    """
    df = pd.read_csv(path_in, sep="\t")

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes dans {path_in} : {missing}")

    df["standardization"] = (df["Reads_Mapped"] * df["gene_length"]) / df["Length"]
    df.to_csv(path_out, sep="\t", index=False)

    return len(df)

# --- Exécution ---
if __name__ == "__main__":

    path_in     = snakemake.input.data
    path_out    = snakemake.output.stand

    n = calculer_score_kegg(path_in, path_out)
    print(f"✓ {n} lignes traitées -> {path_out}")