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
path_out = snakemake.output.taxaname
path_taxonomy = snakemake.input.taxonomy
path_report = snakemake.output.step_checking

def load_taxonomy(path_taxonomy):
    """ Charge le TSV de taxonomie {ID: Lineage} """
    tax_dict = {}
    if not os.path.exists(path_taxonomy):
        print(f"⚠️ Warning: Fichier taxonomie introuvable.")
        return tax_dict
    
    with open(path_taxonomy, 'r', encoding='utf-8') as f:
        for line in f:
            cols = line.strip().split('\t')
            #Format : [0] 'C' | [1] 'ID' | [2] 'TaxID' | [3] 'Lineage'
            if len(cols) >= 4:
                contig_id = cols[1].strip()
                lineage = cols[3].strip()
                # On nettoie le dernier ';' s'il existe
                if lineage.endswith(';'):
                    lineage = lineage[:-1]
                tax_dict[contig_id] = lineage
    return tax_dict

def run_annotation(path_in, path_out, path_taxonomy, tax_dict):
    # 1. Lecture et Annotation
    data_to_sort = []
    tax_dict = load_taxonomy(path_taxonomy)
    
    try:
        with open(path_in, 'r', encoding='utf-8') as f_in:

            for line in f_in:
                cols = line.strip().split('\t')
                if len(cols) < 4: continue
                    
                contig_id = cols[0].strip()
                lineage = tax_dict.get(contig_id, "Unclassified")
                    
                try:
                    rpkm_val = float(cols[3])
                    # On crée la nouvelle ligne enrichie
                    new_line = f"{line.strip()}\t{lineage}\n"
                    data_to_sort.append((rpkm_val, new_line))
                except ValueError: 
                    continue

        # 3. Tri par RPKM (Décroissant)
        data_to_sort.sort(key=lambda x: x[0], reverse=True)
        final_lines = [item[1] for item in data_to_sort]

        # 3. Écriture du fichier final annoté
        with open(path_out, 'w', encoding='utf-8') as f_out:
            f_out.writelines(final_lines)
    
           # Report a 10 lines sample
        with open(path_report, 'w', encoding='utf-8') as f_rep:
            f_rep.write(f"--- APERÇU ALÉATOIRE : {snakemake.wildcards.sample} ---\n")
            f_rep.write("Contig_ID\tLength\tReads\tRPKM\tTaxonomy\n")
            if len(final_lines) > 0:
                sample_size = min(10, len(final_lines))
                random_sample = random.sample(final_lines, sample_size)
                f_rep.writelines(random_sample)
            else:
                f_rep.write("Aucune ligne n'a survécu au filtrage.")

    except Exception as e:
        print(f"❌ Erreur : {e}")
        return None

# --- EXÉCUTION ---
print(f"📊 Analyse de {path_in} ")
if run_annotation(path_in[0], path_in[1], path_out):
    print(f"✅ Annotation réussie -> {os.path.basename(path_out)}")
else:
    # On lève une erreur pour que Snakemake sache que la règle a échoué
    raise RuntimeError(f"Échec de l'annotation pour {os.path.basename(path_in)}")