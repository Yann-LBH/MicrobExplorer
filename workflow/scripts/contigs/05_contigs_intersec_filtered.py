import random
import os

# Shutdown VScode warnings
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(input=SimpleNamespace(), output=SimpleNamespace(), params=SimpleNamespace(), wildcards=SimpleNamespace())
    
# --- RÉCUPÉRATION DES VARIABLES SNAKEMAKE ---
path_in = snakemake.input.data
path_out = snakemake.output.rpkm
path_report = snakemake.output.step_checking

def intersect_filter(path_abund, path_rpkm, path_out):
    """ Garde les lignes du fichier RPKM seulement si le contig est présent dans le fichier Filtered_Contigs. """
    try:
        # 1. Charger les IDs valides du fichier Filtered_Contigs dans un set
        valid_ids = set()
        with open(path_abund, 'r', encoding='utf-8') as f:
            for line in f:
                cols = line.strip().split('\t')
                if cols:
                    valid_ids.add(cols[0])

        kept_lines = []
        # 2. Lire le fichier RPKM et ne garder que ceux présents dans le set
        count_final = 0
        with open(path_rpkm, 'r', encoding='utf-8') as f_in, \
             open(path_out, 'w', encoding='utf-8') as f_out:

            for line in f_in:
                cols = line.strip().split('\t')
                if not cols:
                    continue
                
                contig_id = cols[0]
                
                # On écrit si c'est l'ID est dans les deux fichiers OU si c'est l'en-tête
                if contig_id in valid_ids or contig_id.lower().startswith("contig_id"):
                    f_out.write(line)
                    if not contig_id.lower().startswith("contig_id"):
                        count_final += 1

       # Report a 10 lines sample
        with open(path_report, 'w', encoding='utf-8') as f_rep:
            f_rep.write(f"--- APERÇU ALÉATOIRE : {snakemake.wildcards.sample} ---\n")
            f_rep.write("Contig_ID\tLength\tReads_Mapped\tRPKM\n")
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
print(f"📊 Analyse de {path_in} ")
if intersect_filter(path_in[0], path_in[1], path_out):
    print(f"✅ Intersection réussie -> {os.path.basename(path_out)}")
else:
    # On lève une erreur pour que Snakemake sache que la règle a échoué
    raise RuntimeError(f"Échec de l'intersection pour {os.path.basename(path_in)}")
