################################################################################
# Project : "MicrobExplorer"
# Script: "Merge KEGG with patwahy levels"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging
import pandas as pd

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

FILL_VALUES = {
    "level_1": "Unassigned",
    "level_2": "Unassigned",
    "level_3": "Unassigned",
    "weight": 0,
}


def load_reference(PATHWAY: str) -> pd.DataFrame:
    """Load and prepare the KEGG reference table."""
    ref = pd.read_csv(PATHWAY, sep="\t")
    ref.columns = ref.columns.str.strip()
    ref["ko"] = ref["ko"].astype(str).str.replace("ko:", "", regex=False)
    ref[["level_1", "level_2", "level_3"]] = ref[
        ["level_1", "level_2", "level_3"]
    ].fillna("Unassigned")
    ref["weight"] = ref["weight"].fillna(0)
    return ref


def annotate_with_hierarchy(PATH_IN: str, ref: pd.DataFrame, PATH_OUT: str) -> int:
    """
    KO join -> KEGG hierarchy and weight-adjusted results.
    Returns the number of annotated rows.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    if "kegg" in df.columns:
        df.rename(columns={"kegg": "ko"}, inplace=True)
    df["ko"] = df["ko"].astype(str).str.replace("ko:", "", regex=False)

    merged = df.merge(ref, on="ko", how="left")

    # Ajustement vectorisé sur toutes les colonnes numériques de l'échantillon
    for col in df.select_dtypes(include="number").columns:
        merged[f"adj_{col}"] = merged[col] * merged["weight"]

    merged.fillna(value=FILL_VALUES, inplace=True)
    merged.to_csv(PATH_OUT, sep="\t", index=False)

    return len(merged)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.taxname
    PATHWAY = snakemake.input.pathway

    ref = load_reference(PATHWAY)

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = annotate_with_hierarchy(PATH_IN, ref, PATH_OUT)
    if process:
        logging.info(
            f"[KEGG_MERGE] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_MERGE] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
