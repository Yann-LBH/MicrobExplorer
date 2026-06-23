################################################################################
# Project : "MicrobExplorer"
# Script: "Filter for reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


# All script comments are provided in English as requested.
def filtering_kaiju(PATH_IN, PATH_OUT, COUNT_THRESHOLD):
    """Filters Kaiju output lines based on a minimum read count threshold."""
    try:
        with open(PATH_IN, "r", encoding="utf-8") as f_in, open(
            PATH_OUT, "w", encoding="utf-8"
        ) as f_out:

            # Read and carry over the header line
            header = f_in.readline()
            if header:
                f_out.write(header)

            for line in f_in:
                # Split line by tabs to isolate columns (taxon_id, Reads)
                columns = line.strip().split("\t")

                if len(columns) >= 2:
                    try:
                        # Convert read count column (2nd column, index 1) to integer
                        number_reads = int(columns[1])

                        # Keep line only if it meets or exceeds the threshold
                        if number_reads >= COUNT_THRESHOLD:
                            f_out.write(line)
                    except ValueError:
                        # Skip unexpected text formatting
                        continue

        return True
    except Exception as e:
        print(f"❌ Error processing file {PATH_IN}: {e}")
        return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.filtered
    COUNT_THRESHOLD = snakemake.params.count_threshold

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = filtering_kaiju(PATH_IN, PATH_OUT, COUNT_THRESHOLD)
    if process:
        logging.info(
            f"[READS_FILTER] SUCCESS | Sample: {sample_name} | "
            f"Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_FILTER] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")