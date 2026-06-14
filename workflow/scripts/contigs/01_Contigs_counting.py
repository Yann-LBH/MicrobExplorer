################################################################################
# Project : "MicrobExplorer"
# Script: "Counting contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd


# --- Filtering ---
def filter_contigs(PATH_IN: str, PATH_OUT: str, THRESHOLD: int) -> int:
    """
    Reads the input, filters based on `contig_length` >= `threshold`,
    and writes the output to a TSV file.
    Returns the number of retained contigs.
    """
    df = pd.read_csv(PATH_IN, sep="\t", dtype={"contig_length": int})

    df_filtered = df[df["contig_length"] >= THRESHOLD]

    df_filtered.to_csv(PATH_OUT, sep="\t", index=False)
    return len(df_filtered)


# --- Execution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.raw_data
    PATH_OUT = snakemake.output.counted
    THRESHOLD = int(snakemake.params.length_threshold)

    n_kept = filter_contigs(PATH_IN, PATH_OUT, THRESHOLD)
    print(
        f"✓ CONTIGS : Counting step passed successfully "
        f"-> {n_kept} contigs kept (threshold={THRESHOLD}) "
        f"-> {PATH_OUT}"
    )
