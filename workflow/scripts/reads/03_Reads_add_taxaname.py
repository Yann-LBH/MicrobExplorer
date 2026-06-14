################################################################################
# Project : "MicrobExplorer"
# Script: "Add taxonomy to reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

FINAL_COLS = [
    "tax_id",
    "Reads",
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
    return pd.read_csv(TAXONOMY, sep=";", header=0, quoting=3, dtype={"tax_id": "Int64"})


# ==========================================================================
# 2. Read Kaiju files (format "taxon:reads")
# ==========================================================================
def read_counts(PATH_IN: str) -> pd.DataFrame:
    df = pd.read_csv(PATH_IN, sep=":", header=None, names=["Taxon", "Reads"])
    df["tax_id"] = pd.to_numeric(
        df["Taxon"].str.replace(r"[^0-9]", "", regex=True), errors="coerce"
    ).astype("Int64")
    df["Reads"] = pd.to_numeric(
        df["Reads"].str.replace(r"[^0-9]", "", regex=True), errors="coerce"
    )
    return df[["tax_id", "Reads"]]


# ==========================================================================
# 3. Enrichissement + Top N + Others
# ==========================================================================
def enrich(
    PATH_IN: pd.DataFrame, TAXONOMY: pd.DataFrame, NOISE: list[str], TOP_N: int
) -> pd.DataFrame:

    df = PATH_IN.merge(TAXONOMY, on="tax_id", how="left")
    df = df[~df["domain"].isin(NOISE)]
    df = df[df["scientific_name"].notna()]

    df_top = df.head(TOP_N).copy()
    df_rest = df.iloc[TOP_N:].copy()

    df_others = df_rest.groupby("domain", as_index=False)["Reads"].sum()
    df_others["scientific_name"] = "Others - " + df_others["domain"]
    df_others["tax_id"] = pd.NA

    for col in FINAL_COLS:
        if col not in df_others.columns:
            df_others[col] = "Other"

    return pd.concat([df_top[FINAL_COLS], df_others[FINAL_COLS]], ignore_index=True)


# ==========================================================================
if __name__ == "__main__":
    
    PATH_IN = read_counts(snakemake.input.data)
    PATH_OUT = (snakemake.output.taxaname)

    TAXONOMY = load_mapping(snakemake.input.taxonomy)
    NOISE = list(snakemake.params.noise)
    TOP_N = int(snakemake.params.top_n)

    df_final = enrich(PATH_IN, TAXONOMY, NOISE, TOP_N)
    df_final.to_csv(PATH_OUT, sep="\t", index=False)
    
    print(
        f"✓ READS : Enrichment step passed successfully -> {snakemake.output.taxaname}"
    )