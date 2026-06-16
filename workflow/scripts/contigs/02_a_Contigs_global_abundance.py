import os
import logging
import pandas as pd

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

if __name__ == "__main__":
    # English comments as requested
    INPUT_PATH = snakemake.input.all_data
    OUTPUT_PATH = snakemake.output.global_abundance

    global_counts = pd.Series(dtype=int)

    # Sum abundances across all files
    for f in INPUT_PATH:
        chunk = pd.read_csv(f, sep="\t", usecols=[0, 2])
        contig_col = chunk.columns[0]
        count_col = chunk.columns[1]

        file_counts = chunk.groupby(contig_col)[count_col].sum()
        global_counts = global_counts.add(file_counts, fill_value=0)

    # Save the reference global table
    global_counts.to_csv(
        OUTPUT_PATH, sep="\t", header=["total_abundance"], index_label="contig_id"
    )

    print(f"✓ CONTIGS : Global abundance step passed successfully -> {OUTPUT_PATH}")
