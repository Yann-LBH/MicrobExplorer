################################################################################
# Project : "MicrobExplorer"
# Script: " utils QC Wrapper"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

from utils_qc_Reads import run_full_qc as run_reads
from utils_qc_Contigs import run_full_qc as run_contigs
#from utils_qc_Kegg import run_full_qc as run_kegg

SOURCE = getattr(snakemake.wildcards, "source", None)

if SOURCE == "reads":
    print("--> Launching Quality Control Pipeline for READS")
    run_reads(snakemake)
elif SOURCE == "contigs":
    print("--> Launching Quality Control Pipeline for CONTIGS")
    run_contigs(snakemake)
#elif SOURCE == "kegg":
#    run_kegg()
else:
    raise ValueError(f"Unsupported source wildcard: {SOURCE}")
