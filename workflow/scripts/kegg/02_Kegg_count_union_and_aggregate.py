################################################################################
# Project : "MicrobExplorer"
# Script: "Union between counted contigs and kegg number extracted,
#           Normalization 1/N - Aggregate by kegg number for Deseq"
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


def union_kegg(PATH_IN: str, COUNTS: str, PATH_OUT: str) -> tuple[int, pd.DataFrame]:
    """Union between KEGG annotations and contig counts for a sample.

    Applies 1/N normalization on read counts for genes with multiple KEGG
    annotations. Returns the number of rows in the union and the resulting DataFrame.
    """
    # Read KEGG extraction file from script 1
    df_kegg = pd.read_csv(PATH_IN, sep="\t", keep_default_na=False)
    df_kegg["contig_id"] = df_kegg["contig_id"].astype(str).str.strip()

    # Read contig counts file
    df_counts = pd.read_csv(COUNTS, sep="\t")
    df_counts["contig_id"] = df_counts["contig_id"].astype(str).str.strip()

    if "contig_id" not in df_counts.columns:
        raise KeyError(f"Column 'contig_id' missing in {COUNTS}")

    # Calculate N: number of KO entries associated with each unique gene
    df_kegg["n_ko"] = df_kegg.groupby("contig_id")["contig_id"].transform("count")

    # Inner join between counts and KEGG annotations
    df_out = df_counts.merge(df_kegg, on="contig_id", how="inner")

    if df_out.empty:
        logging.warning(f"No matching contigs found for {COUNTS}")
        df_out.to_csv(PATH_OUT, sep="\t", index=False)
        return 0, df_out

    # Identify numeric count columns from COUNTS to apply 1/N normalization
    metadata_cols = {
        "key",
        "locus_tag",
        "contig_id",
        "kegg_id",
        "gene_length",
        "contig_length",
        "n_ko",
    }
    count_cols = [
        col
        for col in df_out.columns
        if col not in metadata_cols
        and pd.api.types.is_numeric_dtype(df_out[col])
    ]

    # Apply 1/N division on count columns
    for col in count_cols:
        df_out[col] = df_out[col] / df_out["n_ko"]

    # Remove temporary column used for normalization
    df_out.drop(columns=["n_ko"], inplace=True)

    # Export normalized output
    df_out.to_csv(PATH_OUT, sep="\t", index=False)
    return len(df_out), df_out

def aggregate_by_ko(df: pd.DataFrame, MATRIX_DESEQ: str) -> int:
    """Aggregates count values by KEGG KO identifier (kegg_id) for Deseq2."""

    if df.empty or "kegg_id" not in df.columns:
        logging.warning("Empty DataFrame or missing 'kegg_id' column")
        pd.DataFrame(columns=["kegg_id"]).to_csv(
            MATRIX_DESEQ, sep="\t", index=False
        )
        return 0
    
    # Identify numeric count columns dynamically
    count_cols = ["read_mapped"]

    if not count_cols:
        logging.warning("No numeric count columns found for DESeq matrix.")
        return 0

    # Group by KO (including "NA") and sum counts across all count columns
    df_agg = (
        df.groupby("kegg_id", as_index=False)[count_cols]
        .sum()
        .sort_values(by=count_cols[0], ascending=False)
    )

    # Round counts to integers for DESeq2 compatibility
    df_agg[count_cols] = df_agg[count_cols].round(0).astype(int)
    final_cols = ["kegg_id"] + count_cols
    df_agg = df_agg[final_cols]

    df_agg.to_csv(MATRIX_DESEQ, sep="\t", index=False)
    return len(df_agg)

# --- Execution ---
if __name__ == "__main__":
    PATH_IN = snakemake.input.data
    COUNTS = snakemake.input.counted
    PATH_OUT = snakemake.output.union
    MATRIX_DESEQ = snakemake.output.matrix_deseq

    sample_name = getattr(
        snakemake.wildcards, "sample", os.path.basename(PATH_IN)
    )
    # Step 1: Union (returns count and DataFrame in memory)
    n_union, df_union = union_kegg(PATH_IN, COUNTS, PATH_OUT)

    # Step 2: Aggregation using the in-memory DataFrame
    if n_union > 0:
        logging.info(
            f"[KEGG_UNION] SUCCESS | Sample: {sample_name} | Rows: {n_union}"
        )

        n_agg = aggregate_by_ko(df_union, MATRIX_DESEQ)
        logging.info(
            f"[KEGG_MATRIX] SUCCESS | Sample: {sample_name} | Aggregated KOs: {n_agg}"
        )
    else:
        logging.error(
            f"[KEGG_UNION] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )
        raise RuntimeError(f"Filtering failed for {sample_name}")