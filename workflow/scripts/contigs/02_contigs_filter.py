import random
import os

# --- RÉCUPÉRATION DES VARIABLES SNAKEMAKE ---
path_in = snakemake.input.data
# On identifie le fichier spécifique à cet échantillon (wildcard {sample})
current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
path_out = snakemake.output.filtered
path_report = snakemake.output.step_checking
min_total = snakemake.params.reads_threshold

def get_total_abundance(files):
    """ Scanne tous les fichiers pour calculer la somme des reads par contig. """
    global_counts = {}
    
    for path in files:
        with open(path, 'r', encoding='utf-8') as f:           
            for line in f:
                columns = line.strip().split('\t')
                if len(columns) >= 3:
                    try:
                        contig_id = columns[0]
                        reads_mapped = int(columns[2]) #mapped only
                        global_counts[contig_id] = global_counts.get(contig_id, 0) + reads_mapped
                    except ValueError:
                        continue
    return global_counts

def filter_by_global_abundance(current_file, path_out, global_counts, min_total):
    """ Écrit le fichier en vérifiant l'abondance dans le dictionnaire global. """
    kept_lines = []
    try:
        with open(current_file, 'r', encoding='utf-8') as f_in, \
             open(path_out, 'w', encoding='utf-8') as f_out:

            for line in f_in:
                columns = line.strip().split('\t')
                if len(columns) >= 1:
                    contig_id = columns[0]
                    # On garde si c'est l'en-tête OU si l'abondance totale >= 100
                    if contig_id in global_counts:
                        if global_counts[contig_id] >= min_total:
                            f_out.write(line)
                            kept_lines.append(line)

    # Report a 10 lines sample
        with open(path_report, 'w', encoding='utf-8') as f_rep:
            f_rep.write(f"--- APERÇU ALÉATOIRE : {snakemake.wildcards.sample} ---\n")
            f_rep.write(f"Seuil abondance globale : {min_total}\n")
            f_rep.write("Contig_ID\tLength\tReads_Mapped\tReads_Unmapped\n")
            if len(kept_lines) > 0:
                sample_size = min(10, len(kept_lines))
                random_sample = random.sample(kept_lines, sample_size)
                f_rep.writelines(random_sample)
            else:
                f_rep.write("Aucune ligne n'a survécu au filtrage.")

        return True
    except Exception as e:
        print(f"❌ Erreur : {e}")
        return False


# --- EXÉCUTION ---
print(f"📊 Calcul de l'abondance globale sur {len(path_in)} fichiers...")
# Calcul de l'abondance total (nécessaire pour le seuil)
abundance_dict = get_total_abundance(path_in)

print(f"📊 Analyse de {current_file} (Seuil: {min_total})...")
if filter_by_global_abundance(current_file, path_out, abundance_dict, min_total):
    print(f"✅ Filtré avec succès -> {os.path.basename(path_out)}")
else:
    # On lève une erreur pour que Snakemake sache que la règle a échoué
    raise RuntimeError(f"Échec du filtrage pour {current_file}")