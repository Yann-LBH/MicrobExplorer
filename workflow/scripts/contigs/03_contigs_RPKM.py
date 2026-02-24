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

def calculate_rpkm_for_file(path_in, path_out):
    """
    Calcule le RPKM pour chaque ligne d'un fichier.
    Formule : (Reads_Mappés * 10^9) / (Taille_Contig * Total_Reads_Mappés_Echantillon)
    """
    try:
        with open(path_in, 'r', encoding='utf-8') as f:        
            content = f.readlines()
            # 1. Calcul du facteur N (Somme totale des reads mappés dans cet échantillon)
            total_mapped_sample = 0
            valid_data = []
        
            for line in content:
                columns = line.strip().split('\t')
                if len(columns) >= 3:
                    try:
                        mapped = int(columns[2])
                        total_mapped_sample += mapped
                        valid_data.append(columns)
                    except ValueError:
                        continue # Ignore les headers
        
        if total_mapped_sample == 0:
            return False

        # 2. Calcul et écriture du RPKM
        kept_lines = []
        with open(path_out, 'w', encoding='utf-8') as f_out:
            # En-tête
            
            for columns in valid_data:
                contig_id = columns[0]
                length = int(columns[1])
                reads = int(columns[2])
                
                # Application de la formule
                rpkm = (reads * 10**9) / (length * total_mapped_sample)
                
                f_out.write(f"{contig_id}\t{length}\t{reads}\t{rpkm:.4f}\n")
                
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
if calculate_rpkm_for_file(path_in, path_out):
    print(f"✅ Calcul RPKM réussi -> {os.path.basename(path_out)}")
else:
    # On lève une erreur pour que Snakemake sache que la règle a échoué
    raise RuntimeError(f"Échec du calcul RPKM pour {os.path.basename(path_in)}")
