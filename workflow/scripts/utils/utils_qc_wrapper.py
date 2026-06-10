################################################################################
# Project : "MicrobExplorer"
# Script: " utils QC Wrapper"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd
from utils_qc_Reads import run_full_qc as run_reads
from utils_qc_Contigs import run_full_qc as run_contigs
#from utils_qc_Kegg import run_full_qc as run_kegg

# English comments as requested
# Retrieve wildcards and step configurations passed by Snakemake
SOURCE       = snakemake.wildcards.source
STEPS_CONFIG = snakemake.params.steps_config

# Initialize the results variable
results = {}

if SOURCE == "reads":
    print("--> Launching Quality Control Pipeline for READS")
    # ✅ FIX: Store the returned dictionary into the 'results' variable
    results = run_reads(snakemake, STEPS_CONFIG)
    
elif SOURCE == "contigs":
    print("--> Launching Quality Control Pipeline for CONTIGS")
    # ✅ FIX: Store the returned dictionary into the 'results' variable
    results = run_contigs(snakemake, STEPS_CONFIG)

else:
    raise ValueError(f"Unsupported source wildcard: {SOURCE}")

# --- Save Results ---
# Convert the unified results dictionary into a structured Parquet dataframe
df_qc = pd.DataFrame.from_dict(results, orient='index')

# Optional: Ensure index has a proper column name before saving
df_qc.index.name = "Sample"

df_qc.to_parquet(snakemake.output.parquet)
print(f"✓ QC metrics successfully written to {snakemake.output.parquet}")