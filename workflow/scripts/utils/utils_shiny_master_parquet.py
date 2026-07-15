################################################################################
# Project : "MicrobExplorer"
# Script: "Create a master parquet file for the shiny app from multiple parquet files"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging
import polars as pl

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def combine_and_save(file_list, output_path):
    """Helper function to combine a list of parquet files using Polars streaming."""
    # Safety check if a category is empty
    if not file_list:
        logging.warning(f"No input files found for {output_path}. Skipping.")
        return False
        
    logging.info(f"Combining {len(file_list)} files into: {output_path}")
    
    lazy_frames = []
    for file in file_list:
        lf = pl.scan_parquet(file)
        # Add a column with the source file name to identify the sample/test in Shiny
        source_name = os.path.splitext(os.path.basename(file))[0]
        lf = lf.with_columns(pl.lit(source_name).alias("source_file"))
        lazy_frames.append(lf)
        return True
        
    # Stream the concatenated result directly to disk
    combined_query = pl.concat(lazy_frames)
    combined_query.sink_parquet(output_path)
    return False

def main():
    # We iterate over the output names defined in the Snakemake rule.
    # Since inputs and outputs share the exact same keys (pca, heatmap, etc.),
    # we can map them dynamically!
    categories = ["stackedbarplot", "heatmap", "pca", "volcano", "physico", "qc"]
    
    for category in categories:
        # Retrieve the list of input files and the output path from Snakemake
        PATH_IN = getattr(snakemake.input, category)
        PATH_OUT = getattr(snakemake.output, category)
        
        process = combine_and_save(PATH_IN, PATH_OUT)

        if process:
            logging.info(f"[Shiny] SUCCESS | {category.upper()} Master Parquet created : {PATH_OUT}")
        else:
            logging.error(f"[Shiny] FAILED  | {category.upper()} Master Parquet creation failed : {PATH_OUT}")

if __name__ == "__main__":
    main()