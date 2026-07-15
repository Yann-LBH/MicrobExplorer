################################################################################
# Project : "MicrobExplorer"
# Script: " utils QC Wrapper"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd
from utils_qc_Contigs import run_full_qc as run_contigs
from utils_qc_Kegg import run_full_qc as run_kegg
from utils_qc_Reads import run_full_qc as run_reads

# Retrieve wildcards and step configurations passed by Snakemake
SOURCE = snakemake.wildcards.source
STEPS_CONFIG = snakemake.params.steps_config
ACTIVE_MODULES = snakemake.params.active_modules

# Hard-secured guard: Ensure the requested source wildcard is actually declared active
if SOURCE not in ACTIVE_MODULES:
    raise ValueError(
        f"The source pipeline '{SOURCE}' is requested by the workflow but is not marked active "
        f"in the global configuration profile. Active modules are: {ACTIVE_MODULES}"
    )

# Initialize the results variable
results = {}

if SOURCE == "reads":
    print("--> Launching Quality Control Pipeline for READS")
    # Store the returned dictionary into the 'results' variable
    results = run_reads(snakemake, STEPS_CONFIG)

elif SOURCE == "contigs":
    print("--> Launching Quality Control Pipeline for CONTIGS")
    # Store the returned dictionary into the 'results' variable
    results = run_contigs(snakemake, STEPS_CONFIG)

elif SOURCE == "kegg":
    print("--> Launching Quality Control Pipeline for KEGG")
    # Store the returned dictionary into the 'results' variable
    results = run_kegg(snakemake, STEPS_CONFIG)

else:
    raise ValueError(f"Unsupported source wildcard: {SOURCE}")

# --- Save Results ---
# Convert the unified results dictionary into a structured Parquet dataframe
df_qc = pd.DataFrame.from_dict(results, orient="index")

# Ensure index has a proper column name before saving
df_qc.index.name = "Sample"

# Bulletproof Case Standardization: Enforce all tracking metric column names to lowercase
# to prevent any case-mismatch breaks inside downstream R parsing environments.
df_qc.columns = [str(col).lower() for col in df_qc.columns]

# Write DataFrame out to target Parquet file stream
df_qc.to_parquet(snakemake.output.parquet)
print(f"✓ QC metrics successfully written to {snakemake.output.parquet}")
