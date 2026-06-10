################################################################################
# Project : "MicrobExplorer"
# Script: "Standardization by weighted abundance"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################


import pandas as pd

REQUIRED_COLS = {"Reads_Mapped", "gene_length", "Length"}

def calculer_score_kegg(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Calcule l'abondance pondérée : (Reads_Mapped * gene_length) / Length
    Retourne le nombre de lignes traitées.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes dans {PATH_IN} : {missing}")

    df["standardization"] = (df["Reads_Mapped"] * df["gene_length"]) / df["Length"]
    df.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df)

# --- Exécution ---
if __name__ == "__main__":

    PATH_IN     = snakemake.input.data
    PATH_OUT    = snakemake.output.stand

    n = calculer_score_kegg(PATH_IN, PATH_OUT)
    
    print(f"✓ KEGG : Standardization step passed successfully -> {n} lines processed -> {PATH_OUT}")