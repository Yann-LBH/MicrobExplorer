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

def load_reference(taxonomy: str) -> pd.DataFrame:
    """Charge et prépare la table de référence KEGG."""
    ref = pd.read_csv(taxonomy, sep="\t")
    ref.columns = ref.columns.str.strip()
    ref["ko"] = ref["ko"].astype(str).str.replace("ko:", "", regex=False)
    ref[["level_1", "level_2", "level_3"]] = ref[["level_1", "level_2", "level_3"]].fillna("Unassigned")
    ref["weight"] = ref["weight"].fillna(0)
    return ref

def annotate_with_hierarchy(path_in: str, ref: pd.DataFrame, path_out: str) -> int:
    """
    Jointure KO -> hiérarchie KEGG et ajustement pondéré par weight.
    Retourne le nombre de lignes annotées.
    """
    df = pd.read_csv(path_in, sep="\t")

    if "kegg" in df.columns:
        df.rename(columns={"kegg": "ko"}, inplace=True)
    df["ko"] = df["ko"].astype(str).str.replace("ko:", "", regex=False)

    merged = df.merge(ref, on="ko", how="left")

    # Ajustement vectorisé sur toutes les colonnes numériques de l'échantillon
    for col in df.select_dtypes(include="number").columns:
        merged[f"adj_{col}"] = merged[col] * merged["weight"]

    merged.fillna(value=FILL_VALUES, inplace=True)
    merged.to_csv(path_out, sep="\t", index=False)

    return len(merged)

# --- Exécution ---
if __name__ == "__main__":

    path_in     = snakemake.input.data
    path_out    = snakemake.output.taxname
    taxonomy    = snakemake.input.taxonomy

    ref = load_reference(taxonomy)
    n   = annotate_with_hierarchy(path_in, ref, path_out)
    print(f"✓ {n} lignes annotées -> {path_out}")