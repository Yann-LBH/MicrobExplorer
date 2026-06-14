################################################################################
# Project : "MicrobExplorer"
# Script: "Calcul RPKM for contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd


def calculate_rpkm(PATH_IN: str, PATH_OUT: str) -> int:
    """
    Calcul the contigs RPKM
    Formula : (Reads * 10^9) / (Length * Total_Mapped_Reads)
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    required_cols = {"Reads_Mapped", "Length"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes dans {PATH_IN} : {missing}")

    total_mapped = df["Reads_Mapped"].sum()

    if total_mapped == 0:
        print(f"⚠️  Total mapped reads = 0 pour {PATH_IN}, RPKM mis à NaN.")
        df["RPKM"] = float("nan")
        df.to_csv(PATH_OUT, sep="\t", index=False)
        return 0

    df["RPKM"] = ((df["Reads_Mapped"] * 1e9) / (df["Length"] * total_mapped)).round(4)
    df.to_csv(PATH_OUT, sep="\t", index=False)
    return len(df)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.rpkm

    n_rows = calculate_rpkm(PATH_IN, PATH_OUT)

    print(
        f"✓ CONTIGS : RPKM step passed successfully -> {n_rows} contigs processed -> {PATH_OUT}"
    )
