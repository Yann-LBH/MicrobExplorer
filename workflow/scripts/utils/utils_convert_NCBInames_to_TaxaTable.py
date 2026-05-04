import zipfile
import csv
from types import SimpleNamespace

# --- Compatibility stub (dev only) ---
try:
    snakemake
except NameError:
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(zip="data/taxonomy/new_taxdump.zip"),
        output = SimpleNamespace(csv="data/taxonomy/mapping_taxons.csv"),
        params = SimpleNamespace(dmp_name="rankedlineage.dmp")
    )

HEADER = ["tax_id", "scientific_name", "domain", "kingdom",
          "phylum", "class", "order", "family", "genus", "species"]

def extract_taxonomy(zip_path: str, dmp_name: str, out_path: str) -> int:
    """
    Extrait rankedlineage.dmp depuis le ZIP NCBI et écrit un CSV de correspondance.
    Retourne le nombre de lignes écrites.
    """
    n = 0
    with zipfile.ZipFile(zip_path) as z, \
         z.open(dmp_name) as f_in, \
         open(out_path, "w", encoding="utf-8", newline="") as f_out:

        writer = csv.writer(f_out, delimiter=";")
        writer.writerow(HEADER)

        for raw in f_in:
            parts = raw.decode("utf-8").split("|")
            if len(parts) < 10:
                continue

            parts = [p.strip().replace('"', '') for p in parts]
            tax_id, sci_name = parts[0], parts[1]

            # Rangs NCBI ordre inversé (species→kingdom), index 2-9
            ranks = parts[2:10]

            # Contig absent dans un rang -> remplacé par sci_name
            ranks = [r if r else sci_name for r in ranks]

            # Réordonne : domain→kingdom→…→species
            writer.writerow([tax_id, sci_name] + ranks[::-1])
            n += 1

    return n

# --- Exécution ---
if __name__ == "__main__":
    n = extract_taxonomy(
        zip_path = snakemake.input.zip,
        dmp_name = snakemake.params.dmp_name,
        out_path = snakemake.output.csv
    )
    print(f"✓ {n} taxons extraits -> {snakemake.output.csv}")