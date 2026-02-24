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
current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
path_out = snakemake.output.rpkm
path_report = snakemake.output.step_checking
rpkm_threshold = snakemake.params.rpkm_threshold


def get_total_RPKM(path_in):
    """ Scanne tous les fichiers pour verifier le RPKM par contig. """
    min_rpkm_counts = {}
    first_file = True
    
    for path in path_in:
        current_file_contigs = {}
        with open(path, 'r', encoding='utf-8') as f:              
            for line in f:
                colums = line.strip().split('\t')
                if len(colums) >= 4:
                    try:
                        contig_id = colums[0]
                        rpkm = float(colums[3]) # On ne prend que les mappés
                        current_file_contigs[contig_id] = rpkm               
                    except ValueError:
                        continue

        if first_file:
            min_rpkm_counts = current_file_contigs
            first_file = False
        else:
                # On met à jour le dictionnaire minimal
                # 1. On vérifie les contigs communs
                for cid in list(min_rpkm_counts.keys()):
                    if cid in current_file_contigs:
                        min_rpkm_counts[cid] = min(min_rpkm_counts[cid], current_file_contigs[cid])
                    else:
                        # Absent de ce fichier = le minimum tombe à 0
                        min_rpkm_counts[cid] = 0.0

    return min_rpkm_counts

def filter_rpkm_contigs(path_in, path_out, min_rpkm_dict, rpkm_threshold):
    """ Écrit le fichier en vérifiant le RPKM dans le dictionnaire global. """
    kept_lines = []
    try:
        with open(path_in, 'r', encoding='utf-8') as f_in, \
             open(path_out, 'w', encoding='utf-8') as f_out:

            for line in f_in:
                colums = line.strip().split('\t')
                if len(colums) >= 1:
                    contig_id = colums[0]

                    if min_rpkm_dict.get(contig_id, 0.0) >= rpkm_threshold:
                        f_out.write(line)
                        kept_lines.append(line)

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
print(f"📊 RPKM filtering {len(path_in)} fichiers...")

print(f"📊 Analyse de {path_in} ")
if filter_rpkm_contigs(path_in, path_out, min_threshold=0.5):
    print(f"✅ Filtrage RPKM réussi -> {os.path.basename(path_out)}")
else:
    # On lève une erreur pour que Snakemake sache que la règle a échoué
    raise RuntimeError(f"Échec du filtrage RPKM pour {os.path.basename(path_in)}")