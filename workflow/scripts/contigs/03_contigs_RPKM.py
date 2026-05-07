################################################################################
# Project : "MicrobExplorer"
# Script: "Calcul RPKM for contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

def calculate_rpkm(path_in: str, path_out: str) -> int:
    """
    Calcul the contigs RPKM
    Formula : (Reads * 10^9) / (Length * Total_Mapped_Reads)
    """
    df = pd.read_csv(path_in, sep="\t")

    required_cols = {"Reads_Mapped", "Length"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes dans {path_in} : {missing}")

    total_mapped = df["Reads_Mapped"].sum()

    if total_mapped == 0:
        print(f"⚠️  Total mapped reads = 0 pour {path_in}, RPKM mis à NaN.")
        df["RPKM"] = float("nan")
        df.to_csv(path_out, sep="\t", index=False)
        return 0

    df["RPKM"] = ((df["Reads_Mapped"] * 1e9) / (df["Length"] * total_mapped)).round(4)
    df.to_csv(path_out, sep="\t", index=False)
    return len(df)

# --- Exécution ---
if __name__ == "__main__":
    path_in  = snakemake.input.data
    path_out = snakemake.output.rpkm

    n_rows = calculate_rpkm(path_in, path_out)
    print(f"✓ {n_rows} contigs traités -> {path_out}")