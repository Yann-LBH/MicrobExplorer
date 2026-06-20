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


def filtering_kaiju(PATH_IN, PATH_OUT, COUNT_THRESHOLD):
    """Filters Kaiju output lines based on a minimum read count threshold."""
    try:
        with open(PATH_IN, "r", encoding="utf-8") as f_in, open(
            PATH_OUT, "w", encoding="utf-8"
        ) as f_out:

            for line in f_in:
                # Split line by spaces to isolate columns
                columns = line.strip().split(" ")

                if len(columns) >= 4:
                    try:
                        # Convert read count column to integer
                        number_reads = int(columns[3])

                        # Keep line only if it meets or exceeds the threshold
                        if number_reads >= COUNT_THRESHOLD:
                            f_out.write(line)
                    except ValueError:
                        # Skip header lines or unexpected text formatting
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
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_FILTER] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
