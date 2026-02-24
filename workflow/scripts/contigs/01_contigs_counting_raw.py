import random
import os

# --- RÉCUPÉRATION DES VARIABLES SNAKEMAKE ---
path_in = snakemake.input.raw_data
path_out = snakemake.output.counts
min_size = snakemake.params.length_threshold


def filter_tsv_contigs(path_in, path_out, min_size):
    """
    Filtre les fichiers TSV : ne garde que les lignes où la 2ème colonne (taille) >= min_size.
    """
            
    try: 
        os.makedirs(os.path.dirname(path_out), exist_ok=True)

        with open(path_in, 'r', encoding='utf-8') as f_in, \
                open(path_out, 'w', encoding='utf-8') as f_out:
                
            lines = f_in.readlines()
            # On ignore la dernière ligne du fichier TSV
            content = lines[:-1] if len(lines) > 1 else lines
            f_out.write("Contig_ID\tLength\tReads_Mapped\tReads_Unmapped\n")

            for line in content:
                line = line.strip()
                if not line:
                    continue
                
                columns = line.split('\t')
                
                if len(columns) >= 2:
                    try:
                        # Colonne 2 (index 1) : Taille du contig
                        size = int(columns[1])
                        if size >= min_size:
                            f_out.write(line + "\n")
                    except ValueError:
                        # Au cas où la ligne est un en-tête ou du texte
                        continue

        return True
    except Exception as e:
        print(f"❌ Erreur sur le fichier {path_in} : {e}")
        return False

# --- EXÉCUTION ---
print(f"📊 Analyse de {path_in} (Seuil: {min_size})...")
if filter_tsv_contigs(path_in, path_out, min_size):
    print(f"✅ Filtré avec succès -> {path_out}")
else:
    # On lève une erreur pour que Snakemake sache que la règle a échoué
    raise RuntimeError(f"Échec du filtrage pour {path_in}")
