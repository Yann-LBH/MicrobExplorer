################################################################################
# Project : "MicrobExplorer"
# Script: "Create matrix for deseq and Phyloseq for reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os
import logging
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def generate_reads_matrix(
    PATH_IN: str, MATRIX_DESEQ: str, MATRIX_PHYLOSEQ: str
) -> tuple[int, int]:
    """Reads the union contig file and exports raw count and RPKM matrices.

    Returns the number of rows written to each matrix.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    required_cols = {"read_id", "count", "cpm"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns in {PATH_IN}: {missing}")

    if df.empty:
        logging.warning(f"Input file {PATH_IN} is empty.")
        pd.DataFrame(columns=["read_id", "count"]).to_csv(
            MATRIX_DESEQ, sep="\t", index=False
        )
        pd.DataFrame(columns=["read_id", "cpm"]).to_csv(
            MATRIX_PHYLOSEQ, sep="\t", index=False
        )
        return 0, 0

    # 1. Raw Count Matrix for DESeq2
    df_deseq = df[["read_id", "count"]].sort_values(
        by="count", ascending=False
    )
    df_deseq.to_csv(MATRIX_DESEQ, sep="\t", index=False)

    # 2. CPM Matrix for Phyloseq / Abundance
    df_phyloseq = df[["read_id", "cpm"]].sort_values(
        by="cpm", ascending=False
    )
    df_phyloseq.to_csv(MATRIX_PHYLOSEQ, sep="\t", index=False)

    return len(df_deseq), len(df_phyloseq)


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.data
    MATRIX_DESEQ = snakemake.output.matrix_deseq
    MATRIX_PHYLOSEQ = snakemake.output.matrix_phyloseq

    # Report
    sample_name = getattr(
        snakemake.wildcards, "sample", os.path.basename(PATH_IN)
    )
    n_deseq, n_phylo = generate_reads_matrix(
        PATH_IN, MATRIX_DESEQ, MATRIX_PHYLOSEQ
    )

    if n_deseq > 0:
        logging.info(
            f"[READS_MATRIX] SUCCESS | Sample: {sample_name} | "
            f"DESeq2 Rows: {n_deseq} | Phyloseq Rows: {n_phylo}"
        )
    else:
        logging.error(
            f"[READS_MATRIX] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )
        raise RuntimeError(f"Matrix generation failed for {sample_name}")