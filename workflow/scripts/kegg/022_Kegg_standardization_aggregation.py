################################################################################
# Project : "MicrobExplorer"
# Script: "Standardization of aggregation by unique KO per files"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

def agreger_par_ko(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Agrège par code KEGG en sommant la colonne 'standardization'.
    Retourne le nombre de KO uniques.
    """
    df = pd.read_csv(PATH_IN, sep="\t", usecols=["kegg", "standardization"])

    missing = {"kegg", "standardization"} - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes : {missing}")

    df_ko = df.groupby("kegg", as_index=False)["standardization"].sum()
    df_ko.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_ko)

# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.agreg

    n = agreger_par_ko(PATH_IN, PATH_OUT)
    
    print(f"✓ KEGG : Aggregation step passed successfully -> {n} unique KOs -> {PATH_OUT}")