################################################################################
# Project : "MicrobExplorer"
# Script: "Create matrix for deseq and Phyloseq from union contig file"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import logging
import os
import pandas as pd


def generate_contig_matrices(
    PATH_IN: str, MATRIX_DESEQ: str, MATRIX_PHYLOSEQ: str
) -> tuple[int, int]:
    """Reads the union contig file and exports raw count and RPKM matrices.

    Returns the number of rows written to each matrix.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    required_cols = {"contig_id", "raw_count", "rpkm"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Missing required columns in {PATH_IN}: {missing}")

    if df.empty:
        logging.warning(f"Input file {PATH_IN} is empty.")
        pd.DataFrame(columns=["contig_id", "raw_count"]).to_csv(
            MATRIX_DESEQ, sep="\t", index=False
        )
        pd.DataFrame(columns=["contig_id", "rpkm"]).to_csv(
            MATRIX_PHYLOSEQ, sep="\t", index=False
        )
        return 0, 0

    # 1. Raw Count Matrix for DESeq2
    df_deseq = df[["contig_id", "raw_count"]].sort_values(
        by="raw_count", ascending=False
    )
    df_deseq.to_csv(MATRIX_DESEQ, sep="\t", index=False)

    # 2. RPKM Matrix for Phyloseq / Abundance
    df_phyloseq = df[["contig_id", "rpkm"]].sort_values(
        by="rpkm", ascending=False
    )
    df_phyloseq.to_csv(MATRIX_PHYLOSEQ, sep="\t", index=False)

    return len(df_deseq), len(df_phyloseq)


# --- Execution ---
if __name__ == "__main__":
    PATH_IN = snakemake.input.data
    MATRIX_DESEQ = snakemake.output.matrix_deseq
    MATRIX_PHYLOSEQ = snakemake.output.matrix_phyloseq

    sample_name = getattr(
        snakemake.wildcards, "sample", os.path.basename(PATH_IN)
    )

    n_deseq, n_phylo = generate_contig_matrices(
        PATH_IN, MATRIX_DESEQ, MATRIX_PHYLOSEQ
    )

    if n_deseq > 0:
        logging.info(
            f"[CONTIGS_MATRICES] SUCCESS | Sample: {sample_name} | "
            f"DESeq2 Rows: {n_deseq} | Phyloseq Rows: {n_phylo}"
        )
    else:
        logging.error(
            f"[CONTIGS_MATRICES] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )
        raise RuntimeError(f"Matrix generation failed for {sample_name}")