################################################################################
# Project : "MicrobExplorer"
# Script: "Counting contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

# --- Filtering ---
def filter_contigs(path_in: str, path_out: str, threshold: int) -> int:
    """
    Reads the input, filters based on `contig_length` >= `threshold`, and writes the output to a TSV file.
    Returns the number of retained contigs.
    """
    df = pd.read_csv(path_in, sep="\t", dtype={"contig_length": int})

    df_filtered = df[df["contig_length"] >= threshold]

    df_filtered.to_csv(path_out, sep="\t", index=False)
    return len(df_filtered)

# --- Execution ---
if __name__ == "__main__":

    path_in   = snakemake.input.raw_data
    path_out  = snakemake.output.counted
    threshold = int(snakemake.params.length_threshold)

    n_kept = filter_contigs(path_in, path_out, threshold)
    print(f"✓ {n_kept} contigs conservés (seuil={threshold}) -> {path_out}")
