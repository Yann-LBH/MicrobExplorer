################################################################################
# Project : "MicrobExplorer"
# Script: "CPM calculation for reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def calculate_cpm(PATH_IN, PATH_OUT):
    """Calculates Counts Per Million (CPM) for each taxon."""
    try:
        # First pass: Calculate the total number of reads in the file
        total_reads = 0
        data_lines = []
        
        with open(PATH_IN, "r", encoding="utf-8") as f_in:
            header = f_in.readline()  # Skip header
            for line in f_in:
                columns = line.strip().split("\t")
                if len(columns) >= 2:
                    try:
                        reads = int(columns[1])
                        total_reads += reads
                        data_lines.append((columns[0], reads))
                    except ValueError:
                        continue

        # Second pass: Compute CPM and write to output file
        if total_reads > 0:
            with open(PATH_OUT, "w", encoding="utf-8") as f_out:
                # Write new header containing the cpm column
                f_out.write("read_id\tcount\tcpm\n")
                
                for taxon_id, reads in data_lines:
                    # CPM Formula: (reads / total_reads) * 1,000,000
                    cpm = (reads / total_reads) * 1000000
                    f_out.write(f"{taxon_id}\t{reads}\t{cpm:.4f}\n")
            return total_reads
        return False
        
    except Exception as e:
        print(f"❌ Error processing file {PATH_IN}: {e}")
        return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.cpm

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = calculate_cpm(PATH_IN, PATH_OUT)
    if process:
        logging.info(
            f"[READS_CPM] SUCCESS | Sample: {sample_name} | "
            f"Total Reads for CPM: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_CPM] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"CPM calculation failed for {sample_name}")