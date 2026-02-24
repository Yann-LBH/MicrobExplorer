import os

# --- CONFIGURATION ---
FOLDER_FILTERED = "Statistics/2.Filtered_Contigs"
FOLDER_RPKM = "Statistics/4.RPKM_Contigs_Trimmed"
EXIT_FOLDER = "Statistics/52.Final_Intersection_TaxaName"
PATH_TAXONOMY = "A.Taxa_Table_Kaiju/contigs_coassembly-taxNames.tsv" # Ton fichier avec Contig_ID et Lineage

def load_taxonomy(path):
    """ Charge le TSV de taxonomie {ID: Lineage} """
    tax_dict = {}
    if not os.path.exists(path):
        print(f"⚠️ Warning: Fichier taxonomie introuvable.")
        return tax_dict
    
    with open(path, 'r', encoding='utf-8') as f:
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

def intersect_and_annotate(path_abund, path_rpkm, path_out, tax_dict):
    try:
        # 1. IDs valides (Reads cumulés >= 100)
        valid_ids = set()
        with open(path_abund, 'r', encoding='utf-8') as f:
            for line in f:
                cols = line.strip().split('\t')
                if cols: valid_ids.add(cols[0])

        # 2. Lecture, Annotation et Tri
        data_to_sort = []
        header = "Contig_ID\tLength\tReads\tRPKM\tTaxonomy\n"

        with open(path_rpkm, 'r', encoding='utf-8') as f_in:
            # On saute l'ancien header
            next(f_in)
            for line in f_in:
                cols = line.strip().split('\t')
                if len(cols) < 4: continue
                
                contig_id = cols[0]
                if contig_id in valid_ids:
                    # On récupère la taxonomie
                    lineage = tax_dict.get(contig_id, "Unclassified")
                    # On crée la nouvelle ligne enrichie
                    new_line = f"{line.strip()}\t{lineage}\n"
                    
                    try:
                        rpkm_val = float(cols[3])
                        data_to_sort.append((rpkm_val, new_line))
                    except ValueError: continue

        # 3. Tri par RPKM (Décroissant)
        data_to_sort.sort(key=lambda x: x[0], reverse=True)

        # 4. Écriture
        with open(path_out, 'w', encoding='utf-8') as f_out:
            f_out.write(header)
            for _, line in data_to_sort:
                f_out.write(line)
        
        return len(data_to_sort)
    except Exception as e:
        print(f"❌ Erreur : {e}")
        return None

if __name__ == "__main__":
    os.makedirs(EXIT_FOLDER, exist_ok=True)
    tax_dict = load_taxonomy(PATH_TAXONOMY)
    
    for filename in os.listdir(FOLDER_RPKM):
        if filename.endswith(".txt"):
            sample_id = filename.replace("RPKM_Contigs_Trimmed_", "").replace(".txt", "")
            path_rpkm = os.path.join(FOLDER_RPKM, filename)
            path_abund = os.path.join(FOLDER_FILTERED, f"Filtered_Contigs_{sample_id}.txt")
            path_out = os.path.join(EXIT_FOLDER, f"Final_Annotated_{sample_id}.tsv")

            if os.path.exists(path_abund):
                res = intersect_and_annotate(path_abund, path_rpkm, path_out, tax_dict)
                print(f"✅ {sample_id} : {res} contigs annotés et triés.")