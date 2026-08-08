################################################################################
# Project : "MicrobExplorer"
# Script: "Utils Format Megahit Table"
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

# Global definition of standardized lowercase taxonomic ranks
TAX_RANKS = [
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def process_lineage(lineage_str: str) -> list[str]:
    """
    Processes a 7-level lineage string positionally.
    Replaces missing/unwanted levels at rank i with the last known valid rank +
    '_u'.
    """
    UNWANTED = {"unclassified", "unknown", "na", "nan", "null", ""}

    if not isinstance(lineage_str, str) or not lineage_str.strip():
        return ["Unclassified"] * 7

    # Split lineage by semicolon and strip parts
    parts = [p.strip() for p in lineage_str.split(";")]

    # Remove empty trailing strings resulting from a trailing semicolon
    if parts and parts[-1] == "":
        parts.pop()

    ranks = []
    last_valid = "Unclassified"

    # Positional evaluation (0: Domain -> 6: Species)
    for i in range(7):
        if i < len(parts) and parts[i].lower() not in UNWANTED:
            taxon = parts[i]
            ranks.append(taxon)
            last_valid = taxon
        else:
            if last_valid == "Unclassified":
                ranks.append("Unclassified")
            else:
                ranks.append(f"{last_valid}_u")

    return ranks


def load_taxonomy(MEGAHIT: str, TAXONOMY: str) -> pd.DataFrame:
    """Loads Kaiju/Megahit output, formats taxonomy line-by-line, and saves to TSV."""
    # Reading Kaiju format safely (Status, Contig_ID, TaxID, Lineage)
    # Using fillna('') for unclassified contigs ('U') that lack a lineage column
    df = pd.read_csv(
        MEGAHIT,
        sep="\t",
        header=None,
        names=["status", "contig_id", "tax_id", "lineage"],
        usecols=["contig_id", "lineage"],
        dtype=str,
    )
    df["lineage"] = df["lineage"].fillna("")

    # Apply positional processing
    processed_ranks = [process_lineage(l) for l in df["lineage"]]

    expanded_tax = pd.DataFrame(processed_ranks, columns=TAX_RANKS)
    expanded_tax.insert(0, "contig_id", df["contig_id"].values)

    # Save formatted table
    expanded_tax.to_csv(TAXONOMY, sep="\t", index=False)
    return expanded_tax


# --- Exécution ---
if __name__ == "__main__":

    MEGAHIT = str(snakemake.input.megahit)
    TAXONOMY = str(snakemake.output.taxonomy)

    # Report
    process = load_taxonomy(MEGAHIT, TAXONOMY)

    try:
        process = load_taxonomy(MEGAHIT, TAXONOMY)

        if not process.empty:
            logging.info(
                f"[MEGAHIT TABLE CONVERSION] SUCCESS | "
                f"Contigs: {len(process)} | Output: {TAXONOMY}"
            )
        else:
            logging.warning(
                f"[MEGAHIT TABLE CONVERSION] WARNING | File is empty!"
            )

    except Exception as e:
        logging.error(
            f"[MEGAHIT TABLE CONVERSION] FAILED  | Error: {e}"
        )
        raise e