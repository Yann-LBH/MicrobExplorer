################################################################################
# Project : "MicrobExplorer"
# Script: "Data prepared for DEseq2"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

METADATA_COLS = {"contig", "Length", "Contig_ID", "gene_length", "id", "product"}


def sum_per_kegg(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Agrège les counts par KEGG ID.
    Retourne le nombre de KO uniques.
    """
    header = pd.read_csv(PATH_IN, sep="\t", nrows=0).columns
    cols_keep = [c for c in header if c == "kegg" or c not in METADATA_COLS]

    df = pd.read_csv(PATH_IN, sep="\t", usecols=cols_keep)

    if "kegg" not in df.columns:
        raise KeyError(f"Colonne 'kegg' manquante dans {PATH_IN}")

    df_out = df.groupby("kegg").sum(numeric_only=True).reset_index()
    df_out.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_out)


# --- Exécution ---
if __name__ == "__main__":
    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.deseq2

    n = sum_per_kegg(PATH_IN, PATH_OUT)
    print(
        f"✓ KEGG : DESeq2 preparation step passed successfully "
        f"-> {n} unique KOs "
        f"-> {PATH_OUT}"
    )
