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
        PATH_IN, sep="\t", comment="#", header=None, names=gff_cols, dtype=str
    )

    # Filter for specific feature type (e.g., CDS)
    df = df[df["type"] == feature_type].copy()

    # ⚠️ Sécurity Snakemake
    if df.empty:
        logging.warning(
            f"⚠️ No elements of type ‘{feature_type}’ found in {PATH_IN}. "
            f"Creating an empty TSV file for {PATH_OUT}."
        )
        pd.DataFrame(columns=HEADER).to_csv(PATH_OUT, sep="\t", index=False)
        return 0, 0

    if "start" not in df.columns or "end" not in df.columns:
        logging.info(
            f"ℹ️ The ‘start’ or ‘end’ columns are missing in {PATH_IN}. Unable to calculate ‘gene_length’."
        )
        raise KeyError(f"Missing ‘start’/'end' columns in{PATH_IN}")

    df["gene_length"] = (
        pd.to_numeric(df["end"], errors="coerce")
        - pd.to_numeric(df["start"], errors="coerce")
        + 1
    )

    invalid_mask = df["gene_length"].isna() | (df["gene_length"] <= 0)
    if invalid_mask.any():
        n_invalid = int(invalid_mask.sum())
        logging.info(
            f"ℹ️ 'gene_length' missing or invalid for {n_invalid}/{len(df)} genes in {PATH_IN}."
        )
        raise ValueError(
            f"Calcul of 'gene_length' failed for {n_invalid} rows in {PATH_IN}."
        )

    # Conversion explicite en int pour l'export TSV
    df["gene_length"] = df["gene_length"].astype(int)

    # Extract locus_tag (fallback to ID if locus_tag is missing)
    locus_extract = df["attributes"].str.extract(
        r"locus_tag=([^;]+)", expand=False
    )
    id_extract = df["attributes"].str.extract(r"ID=([^;]+)", expand=False)
    df["locus_tag"] = locus_extract.fillna(id_extract).fillna("unknown")
    n_unknown = (df["locus_tag"] == "unknown").sum()
    if n_unknown > 0:
        logging.warning(f"⚠️ {n_unknown}/{len(df)} genes without locus_tag or identifiable ID in {PATH_IN}")

    # Construct global unique gene identifier
    df["key"] = sample_id + "__" + df["locus_tag"]

    # Extract unique KEGG KO identifiers per feature
    df["kegg_id"] = (
        df["attributes"]
        .str.findall(r"K\d{5}")
        .apply(lambda x: sorted(list(set(x))) if x else [])
    )

    total_genes = len(df)
    has_kegg_mask = df["kegg_id"].str.len() > 0
    total_kegg = int(has_kegg_mask.sum())

    # NETTOYAGE MÉMOIRE : Sélection explicite des colonnes de HEADER présentes
    missing = [c for c in HEADER if c not in df.columns and c != "kegg_id"]
    if missing:
        raise ValueError(f"❌ Columns HEADER missing in the parsed GFF : {missing}")
    needed_cols = HEADER 
    df_light = df[needed_cols].copy()
    # Supprime le gros DataFrame d'origine et ses colonnes textuelles lourdes (attributes, etc.)
    del df

    # 5. SÉPARATION : Explode uniquement sur les gènes avec au moins 1 KO
    df_annotated = df_light[has_kegg_mask].explode("kegg_id")

    # Traitement direct des gènes non-annotés (sans explode)
    df_unannotated = df_light[~has_kegg_mask].copy()
    df_unannotated["kegg_id"] = "KO_Unassigned"

    # Concatenation légère
    result_df = pd.concat([df_annotated, df_unannotated], ignore_index=True)

    # Select target columns
    result_df = result_df[HEADER]

    # Export to TSV format
    result_df.to_csv(PATH_OUT, sep="\t", index=False)

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