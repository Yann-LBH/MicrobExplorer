################################################################################
# Project : "MicrobExplorer"
# Script: "Filter for contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

if __name__ == "__main__":
    # English comments as requested
    # Inputs/Outputs are straightforward single files now
    INPUT_PATH = snakemake.input.data
    GLOBAL_ABUNDANCE = snakemake.input.global_table
    OUTPUT_PATH  = snakemake.output.filtered
    THRESHOLD    = int(snakemake.params.abundance_threshold)
    
    # 1. Load the pre-calculated global reference counts
    ref_df = pd.read_csv(GLOBAL_ABUNDANCE, sep="\t", index_col="contig_id")
    
    # 2. Load and filter the single current sample file
    df = pd.read_csv(INPUT_PATH, sep="\t")
    contig_col = df.columns[0]
    
    # Map the target contigs to the pre-computed global values
    df_filtered = df[df[contig_col].map(ref_df["total_abundance"]).fillna(0) >= THRESHOLD]
    
    # 3. Save the individual filtered file
    df_filtered.to_csv(OUTPUT_PATH, sep="\t", index=False)
    
    print(f"✓ CONTIGS : Filtered step passed successfully -> {INPUT_PATH}")