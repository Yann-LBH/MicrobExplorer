import csv
from types import SimpleNamespace

# 1. --- COMPATIBILITY SETTINGS ---
try:
    snakemake
except NameError:
    snakemake = SimpleNamespace(
        input=SimpleNamespace(), 
        output=SimpleNamespace(), 
        params=SimpleNamespace(), 
        wildcards=SimpleNamespace()
    )

# 2. --- JOB DESCRIPTION ---
def filter_tsv_contigs(input_path, output_path, threshold):
    """
    Fonction pure : elle ne dépend d'aucune variable globale.
    Tout ce dont elle a besoin lui est passé en argument.
    """
    try:
        with open(input_path, 'r', encoding='utf-8') as f_in:
            reader = csv.DictReader(f_in, delimiter='\t')
            fieldnames = reader.fieldnames
            kept_rows = [row for row in reader if int(row.get('contig_length', 0)) >= threshold]

        with open(output_path, 'w', encoding='utf-8', newline='') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(kept_rows)
        return True
    
    except Exception as e:
        print(f"❌ Erreur : {e}")
        return False

# 3. --- VARIABLE RETRIEVAL & EXECUTION ---
if __name__ == "__main__":

    path_in = snakemake.input.raw_data
    current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
    path_out = snakemake.output.counts
    min_size = int(snakemake.params.length_threshold)

    print(f"📊 Analyse de {current_file} (Seuil: {min_size})...")

    if filter_tsv_contigs(current_file, path_out, min_size):
        print(f"✅ Filtrage réussi -> {path_out}")
    else:
        raise RuntimeError(f"❌ Échec du traitement pour {current_file}")