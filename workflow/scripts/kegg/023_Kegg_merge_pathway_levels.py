################################################################################
# Project : "MicrobExplorer"
# Script: "Merge KEGG with patwahy levels"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

FILL_VALUES = {
    "level_1": "Unassigned",
    "level_2": "Unassigned",
    "level_3": "Unassigned",
    "weight":  0
}

def load_reference(PATHWAY: str) -> pd.DataFrame:
    """Charge et prépare la table de référence KEGG."""
    ref = pd.read_csv(PATHWAY, sep="\t")
    ref.columns = ref.columns.str.strip()
    ref["ko"] = ref["ko"].astype(str).str.replace("ko:", "", regex=False)
    ref[["level_1", "level_2", "level_3"]] = ref[["level_1", "level_2", "level_3"]].fillna("Unassigned")
    ref["weight"] = ref["weight"].fillna(0)
    return ref

def annotate_with_hierarchy(PATH_IN: str, ref: pd.DataFrame, PATH_OUT: str) -> int:
    """
    Jointure KO -> hiérarchie KEGG et ajustement pondéré par weight.
    Retourne le nombre de lignes annotées.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    if "kegg" in df.columns:
        df.rename(columns={"kegg": "ko"}, inplace=True)
    df["ko"] = df["ko"].astype(str).str.replace("ko:", "", regex=False)

    merged = df.merge(ref, on="ko", how="left")

    # Ajustement vectorisé sur toutes les colonnes numériques de l'échantillon
    for col in df.select_dtypes(include="number").columns:
        merged[f"adj_{col}"] = merged[col] * merged["weight"]

    merged.fillna(value=FILL_VALUES, inplace=True)
    merged.to_csv(PATH_OUT, sep="\t", index=False)

    return len(merged)

# --- Exécution ---
if __name__ == "__main__":

    PATH_IN     = snakemake.input.data
    PATH_OUT    = snakemake.output.taxname
    PATHWAY     = snakemake.input.pathway

    ref = load_reference(PATHWAY)
    n   = annotate_with_hierarchy(PATH_IN, ref, PATH_OUT)
    
    print(f"✓ KEGG : Merge step passed successfully -> {n} lines annotated -> {PATH_OUT}")