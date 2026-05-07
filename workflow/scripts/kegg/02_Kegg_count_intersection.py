################################################################################
# Project : "MicrobExplorer"
# Script: "Intersection between counted contigs and kegg number extracted"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

def intersection_kegg(path_in: str, counts: str, path_out: str) -> int:
    """
    Jointure inner entre les annotations KEGG et les counts d'un échantillon.
    Retourne le nombre de lignes dans l'intersection.
    """
    df_kegg = pd.read_csv(path_in, sep="\t")
    df_kegg["contig"] = df_kegg["contig"].astype(str).str.strip()

    df_counts = pd.read_csv(counts, sep="\t")
    df_counts["Contig_ID"] = df_counts["Contig_ID"].astype(str).str.strip()
    df_counts.drop(columns=["Reads_Unmapped"], errors="ignore", inplace=True)

    if "Contig_ID" not in df_counts.columns:
        raise KeyError(f"Colonne 'Contig_ID' manquante dans {counts}")

    df_out = df_counts.merge(df_kegg, left_on="Contig_ID", right_on="contig", how="inner")

    if df_out.empty:
        print(f"⚠️ Aucun contig commun pour {counts}")

    df_out.to_csv(path_out, sep="\t", index=False)
    return len(df_out)

# --- Exécution ---
if __name__ == "__main__":

    path_in   = snakemake.input.data
    counts = snakemake.input.counted
    path_out = snakemake.output.intersec
 
    n = intersection_kegg(path_in, counts, path_out)
    print(f"✓ {n} contigs en intersection -> {path_out}")