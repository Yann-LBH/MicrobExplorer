################################################################################
# Project : "MicrobExplorer"
# Script: "Standardization of aggregation by unique KO per files"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

def agreger_par_ko(path_in: str, path_out: str) -> int:
    """
    Agrège par code KEGG en sommant la colonne 'standardization'.
    Retourne le nombre de KO uniques.
    """
    df = pd.read_csv(path_in, sep="\t", usecols=["kegg", "standardization"])

    missing = {"kegg", "standardization"} - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes : {missing}")

    df_ko = df.groupby("kegg", as_index=False)["standardization"].sum()
    df_ko.to_csv(path_out, sep="\t", index=False)

    return len(df_ko)

# --- Exécution ---
if __name__ == "__main__":

    path_in = snakemake.input.data
    path_out = snakemake.input.agreg

    n = agreger_par_ko(path_in, path_out)
    print(f"✓ {n} KO uniques -> {path_out}")