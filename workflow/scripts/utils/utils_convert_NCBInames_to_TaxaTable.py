import os
import urllib.request
import logging
import zipfile
import csv

# Configure logging to display time, level, and message properly
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

HEADER = [
    "tax_id",
    "scientific_name",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

# ==========================================================================
# Step 1 : downloading zipfile from NCBI
# ==========================================================================
def download_ncbi_taxonomy(URL: str, ZIP_PATH: str) -> None:
    """
    Downloads the NCBI taxonomy zip file if it does not already exist.
    """
    if not os.path.exists(ZIP_PATH):
        logging.info(f"Downloading file from {URL}...")
        urllib.request.urlretrieve(URL, ZIP_PATH)

# ==========================================================================
# Step 2 : Extraction
# ==========================================================================
def extract_taxonomy(ZIP_PATH: str, DMP_NAME: str, TAXONOMY: str) -> int:
    """
    Extracts taxonomy information from the specified .dmp file within the zip archive 
    and writes it to a TSV file compatible with the pipeline.
    """
    n = 0
    with zipfile.ZipFile(ZIP_PATH) as z, z.open(os.path.basename(DMP_NAME)) as f_in, open(
        TAXONOMY, "w", encoding="utf-8", newline=""
    ) as f_out:

        # Utilisation d'un séparateur TABULATION (\t) pour le format TSV
        writer = csv.writer(f_out, delimiter="\t")
        writer.writerow(HEADER)

        for raw in f_in:
            line = raw.decode("utf-8")
            parts = [p.strip() for p in line.split("|")]

            # Security
            #if len(parts) < 9:
            #    continue

            tax_id = parts[0]
            sci_name = parts[1].replace('"', "")

            # Récupération des rangs de l'index 2 à 10 exclus (comme votre ancien script)
            ranks = [p.replace('"', "") for p in parts[2:10]]
            
            # Inversion directe pour passer de l'ordre NCBI (species->domain) à (domain->species)
            ranks = ranks[::-1]

            # --- Logique de remplissage de l'ancien script ---
            first_empty_found = False
            for i in range(len(ranks)):
                if ranks[i] == "":
                    if not first_empty_found:
                        ranks[i] = sci_name
                        first_empty_found = True
                    else:
                        ranks[i] = f"{sci_name}_u"

            # Écriture de la ligne finale dans le fichier TSV
            writer.writerow([tax_id, sci_name] + ranks)
            n += 1

    return n


# --- Exécution ---
if __name__ == "__main__":

    URL = snakemake.params.url
    DMP_NAME = snakemake.params.dmp_name
    ZIP_PATH = snakemake.output.zip
    TAXONOMY = snakemake.output.taxonomy  

    try:
        logging.info(f"Downloading NCBI taxonomy from {URL}...")
        download_ncbi_taxonomy(URL, ZIP_PATH)
        logging.info(f"Download complete: {ZIP_PATH}")

        logging.info(f"Starting taxonomy extraction from {ZIP_PATH}...")
        lines_written = extract_taxonomy(ZIP_PATH, DMP_NAME, TAXONOMY)
        
        if lines_written > 0:
            logging.info(
                f"[TAXONOMY EXTRACTION] SUCCESS | {lines_written} rows written to {TAXONOMY}"
            )
        else:
            logging.error("[TAXONOMY EXTRACTION] FAILED | No rows processed.")
            raise RuntimeError("Extraction produced an empty file.")

    except Exception as e:
        logging.error(f"[TAXONOMY EXTRACTION] FAILED with error: {str(e)}")
        raise e