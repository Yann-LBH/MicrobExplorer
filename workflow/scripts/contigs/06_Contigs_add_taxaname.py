################################################################################
# Project : "MicrobExplorer"
# Script: "Add taxonomy to contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

def load_taxonomy(TAXONOMY: str) -> pd.Series:
    """
    Loads the taxonomy file.
    Expected format: [0] Class | [1] ID | [2] TaxID | [3] Lineage
    Returns a pandas Series indexed by Contig_ID.
    """
    df = pd.read_csv(TAXONOMY, sep="\t", header=None,
                     usecols=[1, 3], names=["Contig_ID", "Lineage"])
    df["Lineage"] = df["Lineage"].str.rstrip(";")
    return df.set_index("Contig_ID")["Lineage"]

def run_annotation(PATH_IN: str, PATH_OUT: str, taxonomy: pd.Series) -> int:
    """
    Annotates RPKM data with taxonomy, sorts by RPKM descending, and writes the output.
    Returns the number of annotated contigs.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    df["Taxonomy"] = df["Contig_ID"].map(taxonomy).fillna("Unclassified")

    df = df.sort_values("RPKM", ascending=False)
    df.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df)

# --- Exécution ---
if __name__ == "__main__":
    PATH_IN     = snakemake.input.data
    PATH_OUT    = snakemake.output.taxaname
    TAXONOMY    = snakemake.input.taxonomy

    taxonomy = load_taxonomy(TAXONOMY)

    n_kept = run_annotation(PATH_IN, PATH_OUT, taxonomy)
    
    print(f"✓ CONTIGS : Taxonomy annotation step passed successfully -> {n_kept} contigs annotated -> {PATH_OUT}")