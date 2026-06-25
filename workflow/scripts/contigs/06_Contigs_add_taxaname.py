################################################################################
# Project : "MicrobExplorer"
# Script: "Add taxonomy to contigs"
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


# All script comments are provided in English as requested.
import logging
import os
import pandas as pd

# Global definition of standardized lowercase taxonomic ranks
TAX_RANKS = [
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def load_taxonomy(TAXONOMY: str) -> pd.DataFrame:
    """Loads the taxonomy file and expands the semi-colon separated Lineage

    into 8 standardized lowercase columns.
    Expected input format: [0] Class | [1] ID | [2] TaxID | [3] Lineage
    """
    # Load columns 1 (ID) and 3 (Lineage)
    df = pd.read_csv(
        TAXONOMY,
        sep="\t",
        header=None,
        usecols=[1, 3],
        names=["contig_id", "lineage"],
    )

    # Clean trailing semi-colons and split into separate columns
    df["lineage"] = df["lineage"].str.rstrip(";")
    expanded_tax = (
        df["lineage"].str.split(r";\s*", expand=True).fillna("unclassified")
    )

    # Slice or pad columns to match exactly the 8 standard ranks
    num_cols = expanded_tax.shape[1]
    if num_cols < 8:
        for i in range(num_cols, 8):
            expanded_tax[i] = "unclassified"
    elif num_cols > 8:
        expanded_tax = expanded_tax.iloc[:, :8]

    # Assign lowercase taxonomic names as column headers
    expanded_tax.columns = TAX_RANKS
    expanded_tax["contig_id"] = df["contig_id"]

    return expanded_tax.set_index("contig_id")


def run_annotation(PATH_IN: str, PATH_OUT: str, df_tax: pd.DataFrame) -> int:
    """Merges the lowercase expanded taxonomy columns into the incoming RPKM

    dataset.
    """
    df = pd.read_csv(PATH_IN, sep="\t")

    # Ensure columns from step 5 are also strictly lowercase
    df.columns = df.columns.str.lower()

    # Join the expanded taxonomy data based on contig_id
    df_annotated = df.merge(df_tax, on="contig_id", how="left")

    # Fill any missing unmatched contigs with 'unclassified' for all ranks
    df_annotated[TAX_RANKS] = df_annotated[TAX_RANKS].fillna("unclassified")

    # Replace empty strings or whitespace-only strings with 'unclassified'
    for col in TAX_RANKS:
        df_annotated[col] = df_annotated[col].replace(
            r"^\s*$", "unclassified", regex=True
        )

    df_annotated.to_csv(PATH_OUT, sep="\t", index=False)

    return len(df_annotated)


# --- Exécution ---
if __name__ == "__main__":
    PATH_IN = snakemake.input.data
    PATH_OUT = snakemake.output.taxaname
    TAXONOMY = snakemake.input.taxonomy

    taxonomy = load_taxonomy(TAXONOMY)

    # Report
    sample_name = getattr(snakemake.wildcards, "sample", os.path.basename(PATH_IN))
    process = run_annotation(PATH_IN, PATH_OUT, taxonomy)
    if process:
        logging.info(
            f"[CONTIGS_ADD_TAXANAME] SUCCESS | Sample: {sample_name} | "
            ""
            f"Count: {process} | Output: {PATH_OUT}"
        )
    else:
        logging.error(
            f"[CONTIGS_ADD_TAXANAME] FAILED  | Sample: {sample_name} | Input: {PATH_IN}"
        )

        raise RuntimeError(f"Filtering failed for {sample_name}")
