################################################################################
# Project : "MicrobExplorer"
# Script: "Calcul RPKM for Kegg - Aggregate by Kegg number for Phyloseq"
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


def calculate_rpkm(PATH_IN: str, PATH_OUT: str) -> tuple[int, pd.DataFrame]:
    """
    Calcul the kegg RPKM
    Formula : (Read_mapped (contigs) * 10^9) / (Contig_length * Total_Mapped_Reads)
    """
    df_rpkm = pd.read_csv(PATH_IN, sep="\t", keep_default_na=False)

    required_cols = {"contig_length", "read_mapped"}
    missing = required_cols - set(df_rpkm.columns)
    if missing:
        raise KeyError(f"Missing column in {PATH_IN} : {missing}")

    total_mapped = df_rpkm["read_mapped"].sum()

    if total_mapped == 0:
        print(f"⚠️  Total mapped reads = 0 for {PATH_IN}, RPKM set to NaN.")
        df_rpkm["rpkm"] = float("nan")
        df_rpkm.to_csv(PATH_OUT, sep="\t", index=False)
        return 0, df_rpkm

    df_rpkm["rpkm"] = (
    (df_rpkm["read_mapped"] * 1e9) / (df_rpkm["contig_length"] * total_mapped)
    ).round(4)
    df_rpkm.to_csv(PATH_OUT, sep="\t", index=False)
    return len(df_rpkm), df_rpkm

def aggregate_by_ko(df: pd.DataFrame, MATRIX_PHYLOSEQ: str) -> int:
    """Aggregates count values by KEGG KO identifier (kegg_id) for Phyloseq."""

    if df.empty or "kegg_id" not in df.columns:
        logging.warning("Empty DataFrame or missing 'kegg_id' column")
        pd.DataFrame(columns=["kegg_id"]).to_csv(
            MATRIX_PHYLOSEQ, sep="\t", index=False
        )
        return 0
    
    # Identify numeric count columns dynamically
    count_cols = ["rpkm"]
    # Group by KO (including "NA") and sum counts across all count columns
    df_agg = (
        df.groupby("kegg_id", as_index=False)[count_cols]
        .sum()
        .sort_values(by=count_cols[0], ascending=False)
    )

    df_agg.to_csv(MATRIX_PHYLOSEQ, sep="\t", index=False)
    return len(df_agg)


# --- Exécution ---
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.rpkm
    MATRIX_PHYLOSEQ = snakemake.output.matrix_phyloseq

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    # Step 1: Intersection (returns count and DataFrame in memory)
    n_rpkm, df_rpkm = calculate_rpkm(PATH_IN, PATH_OUT)

    # Step 2: Aggregation using the in-memory DataFrame
    if n_rpkm > 0:
        logging.info(
            f"[KEGG_RPKM] SUCCESS | Sample: {sample_name} | Rows: {n_rpkm}"
        )

        n_agg = aggregate_by_ko(df_rpkm, MATRIX_PHYLOSEQ)
        logging.info(
            f"[KEGG_MATRIX] SUCCESS | Sample: {sample_name} | Aggregated KOs: {n_agg}"
        )
    else:
        logging.error(
            f"[KEGG_RPKM] FAILED | Sample: {sample_name} | Input: {PATH_IN}"
        )
        raise RuntimeError(f"Processing failed for {sample_name}")