# ==============================================================================
# PROJECT   : MicrobExplorer
# SCRIPT    : utils_pathway_level_extraction.py
# PURPOSE   : Download KEGG Pathway and create a pathwaytable
# AUTHOR    : Yann Le Bihan
# DATE      : 2026-09-03
# LINK      : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

import logging
import re
import pandas as pd
import requests

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def get_kegg_hierarchy(url: str, output_path: str) -> pd.DataFrame:
    """
    Fetches and parses the KEGG hierarchy (ko00001) to extract
    functional levels and EC numbers, then saves the formatted table to TSV.
    """
    # Pre-compile regex patterns for significant speedup in loops
    re_clean_ab = re.compile(r"^[A-B\s]*\d+\s+")
    re_pathway = re.compile(r"^\d+\s+(.*?)(?:\s*\[(?:PATH|BR):ko\d+\])?$")
    re_ec = re.compile(r"\[EC:(.*?)\]")

    logging.info(f"Fetching data from KEGG API: {url}")
    response = requests.get(url, timeout=30)
    response.raise_for_status()

    hierarchy_data = []
    l1, l2, l3 = "", "", ""

    # Split text into lines once to optimize iteration performance
    lines = response.text.splitlines()

    for line in lines:
        if not line:
            continue

        prefix = line[0]

        if prefix == "A":
            l1 = re_clean_ab.sub("", line).strip()
        elif prefix == "B":
            l2 = re_clean_ab.sub("", line).strip()
        elif prefix == "C":
            # Fast extraction of pathway name
            match_c = re_pathway.search(line[1:].strip())
            l3 = match_c.group(1).strip() if match_c else line[3:].strip()
        elif prefix == "D":
            line_content = line[4:].strip()

            # Extract EC numbers if they exist
            ec_match = re_ec.search(line_content)
            ec_val = ec_match.group(1) if ec_match else ""

            # Clean content for ID and Description extraction
            clean_content = re_ec.sub("", line_content).strip()

            # Efficient split to capture KO identifier and its associated description
            parts = clean_content.split(None, 1)
            ko_id = parts[0]
            ko_desc = parts[1] if len(parts) > 1 else ""

            hierarchy_data.append(
                {
                    "ko": ko_id,
                    "ec_number": ec_val,
                    "level_1": l1,
                    "level_2": l2,
                    "level_3": l3,
                    "gene_description": ko_desc,
                }
            )

    df = pd.DataFrame(hierarchy_data)

    df.to_csv(output_path, sep="\t", index=False)
    return df


# --- Snakemake Execution ---
if __name__ == "__main__":

    URL = str(snakemake.params.url)
    PATHWAY = str(snakemake.output.pathway)

    try:
        process = get_kegg_hierarchy(URL, PATHWAY)

        if not process.empty:
            logging.info(
                f"[PATHWAY EXTRACTION] SUCCESS | Count: {len(process)} entries | Output: {PATHWAY}"
            )
        else:
            logging.warning(
                f"[PATHWAY EXTRACTION] WARNING | Output file is empty: {PATHWAY}"
            )

    except Exception as e:
        logging.error(
            f"[PATHWAY EXTRACTION] FAILED  | URL: {URL} | Error: {e}"
        )
        raise e