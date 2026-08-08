################################################################################
# Project : "MicrobExplorer"
# Script: "Extraction of lines with and without a KEGG number"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import logging
import os
import pandas as pd

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

HEADER = ["key", "locus_tag", "contig_id", "kegg_id", "gene_length"]


def process_gff_kegg(
    PATH_IN: str, PATH_OUT: str, sample_id: str, feature_type: str = "CDS"
) -> tuple[int, int]:
    """Parses a GFF3 file using pandas, extracts global gene IDs and KEGG annotations.

    Duplicates rows for genes with multiple KEGG terms (1 row per KO)
    and keeps genes without KEGG annotations as 'NA'.
    No read counts division is performed at this stage.
    """
    gff_cols = [
        "contig_id",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    # Read GFF3 file skipping comment lines
    df = pd.read_csv(
        PATH_IN, sep="\t", comment="#", header=None, names=gff_cols
    )

    # Filter for specific feature type (e.g., CDS)
    df = df[df["type"] == feature_type].copy()

    if df.empty:
        return 0, 0

    # Calculate gene length
    df["gene_length"] = df["end"] - df["start"] + 1

    # Extract locus_tag (fallback to ID if locus_tag is missing)
    locus_extract = df["attributes"].str.extract(
        r"locus_tag=([^;]+)", expand=False
    )
    id_extract = df["attributes"].str.extract(r"ID=([^;]+)", expand=False)
    df["locus_tag"] = locus_extract.fillna(id_extract).fillna("unknown")

    # Construct global unique gene identifier
    df["key"] = sample_id + "__" + df["locus_tag"]

    # Extract unique KEGG KO identifiers per feature
    df["kegg_id"] = (
        df["attributes"]
        .str.findall(r"K\d{5}")
        .apply(lambda x: sorted(list(set(x))) if x else ["NA"])
    )

    # Explode the KEGG list so each KO gets its own row
    df_exploded = df.explode("kegg_id")

    # Select target columns
    result_df = df_exploded[HEADER]

    # Export to TSV format
    result_df.to_csv(PATH_OUT, sep="\t", index=False)

    total_genes = len(df)
    total_kegg = int((df_exploded["kegg_id"] != "NA").sum())

    return total_genes, total_kegg


# --- Execution ---
if __name__ == "__main__":
    PATH_IN = snakemake.input.raw_data
    PATH_OUT = snakemake.output.extracted

    # Retrieve sample name from Snakemake wildcards or fallback to filename
    sample_name = getattr(
        snakemake.wildcards, "sample", os.path.basename(PATH_IN).split(".")[0]
    )

    total_genes, total_kegg = process_gff_kegg(
        PATH_IN, PATH_OUT, sample_id=sample_name
    )

    if total_genes > 0:
        logging.info(
            f"[KEGG_EXTRACT] SUCCESS | Sample: {sample_name} | "
            f"Total unique CDS: {total_genes} | KEGG terms extracted: {total_kegg} | "
            f"Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[KEGG_EXTRACT] FAILED  | Sample: {sample_name} | "
            f"No CDS features found in Input: {PATH_IN}"
        )
        raise RuntimeError(
            f"Extraction failed: no CDS found for sample {sample_name}"
        )