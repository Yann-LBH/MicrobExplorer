################################################################################
# Project : "MicrobExplorer"
# Script: "Add taxonomy to reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import logging
import pandas as pd

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

FINAL_COLS = [
    "Read_ID",
    "Count",
    "CPM",
    "scientific_name",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


# ==========================================================================
# 1. Loading mapping NCBI
# ==========================================================================
def load_mapping(TAXONOMY: str) -> pd.DataFrame:
    # Directly reading taxonomy assuming the column is already "Read_ID"
    return pd.read_csv(
        TAXONOMY, sep="\t", header=0, quoting=3, dtype={"Read_ID": str}
    )


# ==========================================================================
# 2. Read Kaiju files (TSV format: Read_ID, Count, CPM)
# ==========================================================================
def read_counts(PATH_IN: str) -> pd.DataFrame:
    # Directly loading file without any renaming steps
    return pd.read_csv(
        PATH_IN,
        sep="\t",
        header=0,
        dtype={"Read_ID": str, "Count": int, "CPM": float},
    )


# ==========================================================================
# 3. Enrichment + Top N + Others
# ==========================================================================
def enrich(
    PATH_IN: pd.DataFrame, TAXONOMY: pd.DataFrame, NOISE: list[str], TOP_N: int
) -> pd.DataFrame:

    df = PATH_IN.merge(TAXONOMY, left_on="Read_ID", right_on="tax_id", how="left")
    df = df[~df["domain"].isin(NOISE)]
    df = df[df["scientific_name"].notna()]

    df_top = df.head(TOP_N).copy()
    df_rest = df.iloc[TOP_N:].copy()

    # Summing up both Count and CPM for the remaining entries
    df_others = (
        df_rest.groupby("domain", as_index=False)[["Count", "CPM"]].sum()
    )
    df_others["scientific_name"] = "Others - " + df_others["domain"]
    df_others["Read_ID"] = pd.NA

    for col in FINAL_COLS:
        if col not in df_others.columns:
            df_others[col] = "Other"

    return pd.concat(
        [df_top[FINAL_COLS], df_others[FINAL_COLS]], ignore_index=True
    )


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = read_counts(snakemake.input.data)
    PATH_OUT = snakemake.output.taxaname
    TAXONOMY = load_mapping(snakemake.input.taxonomy)
    NOISE = list(snakemake.params.noise)
    TOP_N = int(snakemake.params.top_n)

    # Report
    sample_name = snakemake.wildcards.sample

    if "Read_ID" in PATH_IN.columns:
        PATH_IN["Read_ID"] = PATH_IN["Read_ID"].astype(str).str.strip()
    if "Read_ID" in TAXONOMY.columns:
        TAXONOMY["Read_ID"] = TAXONOMY["Read_ID"].astype(str).str.strip()

    process = enrich(PATH_IN, TAXONOMY, NOISE, TOP_N)
    process.to_csv(PATH_OUT, sep="\t", index=False)
    if not process.empty:
        logging.info(
            f"[READS_ADD_TAXANAME] SUCCESS | Sample: {sample_name} | "
            f"Count: {len(process)} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[READS_ADD_TAXANAME] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")