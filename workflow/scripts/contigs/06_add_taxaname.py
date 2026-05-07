################################################################################
# Project : "MicrobExplorer"
# Script: "Add taxonomy to contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

def load_taxonomy(taxonomy: str) -> pd.Series:
    """
    Charge le fichier taxonomie.
    Format attendu : [0] Class | [1] ID | [2] TaxID | [3] Lineage
    Retourne une Series indexée par Contig_ID.
    """
    df = pd.read_csv(taxonomy, sep="\t", header=None,
                     usecols=[1, 3], names=["Contig_ID", "Lineage"])
    df["Lineage"] = df["Lineage"].str.rstrip(";")
    return df.set_index("Contig_ID")["Lineage"]

def run_annotation(path_in: str, path_out: str,
                   taxonomy: pd.Series) -> int:
    """
    Annote les données RPKM avec la taxonomie,
    trie par RPKM décroissant et écrit le résultat.
    Retourne le nombre de contigs annotés.
    """
    df = pd.read_csv(path_in, sep="\t")

    df["Taxonomy"] = df["Contig_ID"].map(taxonomy).fillna("Unclassified")

    df = df.sort_values("RPKM", ascending=False)
    df.to_csv(path_out, sep="\t", index=False)

    return len(df)

# --- Exécution ---
if __name__ == "__main__":
    path_in  = snakemake.input.data
    path_out = snakemake.output.taxaname
    taxonomy    = snakemake.output.taxonomy

    taxonomy = load_taxonomy(snakemake.input.taxonomy)

    n_kept = run_annotation(path_in, path_out, taxonomy)
    print(f"✓ {n_kept} contigs annotés -> {path_out}")