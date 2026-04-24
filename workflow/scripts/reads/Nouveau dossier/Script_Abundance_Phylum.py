import os
import zipfile
import shutil

# --- CONFIGURATION (En haut pour être facilement modifiable) ---
ZIP_PATH = "A.Taxa_Table_NCBI/new_taxdump.zip"
NAME_FILE_IN_ZIP = "rankedlineage.dmp"
ROOT_FOLDER = "Data"
EXIT_FOLDER = "Abundance/Phylum"

def phylum_abundance(ZIP_PATH, NAME_FILE_IN_ZIP, ROOT_FOLDER, EXIT_FOLDER):
    # 1. Create output directory if it doesn't exist    
    output_dir = os.path.dirname(EXIT_FOLDER)

    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Dossier créé : {output_dir}")

    print(f"Détermination du phylum des reads depuis {NAME_FILE_IN_ZIP}...")
    tax_to_phylum = {}
    
    with zipfile.ZipFile(ZIP_PATH, 'r') as z:
        with z.open(NAME_FILE_IN_ZIP) as f_in:
            
            for binary_line in f_in:
                parts = [p.strip() for p in binary_line.decode('utf-8').split('|')]

                if len(parts) >= 9:
                    tax_to_phylum[parts[0]] = parts[7] if parts[7] else "Unknown Phylum" #Associe TaxaId au Domain 0 = tax_id, 9 = domain
    
    # 2. Comptage des reads dans le fichier Kaiju
    counter = {}
    total_reads = 0

    with open(ROOT_FOLDER, 'r', encoding='utf-8') as f:
        for line in f:
            total_reads += 1
            column = line.strip().split('\t')
                
            # Kaiju : C = Classified, U = Unclassified
            status = column[0]
                
            if status == "C":
                tax_id = column[2]
                phylum = tax_to_phylum.get(tax_id, "Unknown/Other") #Si le tax_id n'est pas dans le mapping
                counter[phylum] = counter.get(phylum, 0) + 1
            else:
                counter["Unclassified"] = counter.get("Unclassified", 0) + 1
            
    if total_reads > 0:
        sorted_count = sorted(counter.items(), key=lambda x: x[1], reverse=True)
        #Garde les 19 premières lignes seulement ( les plus abondants) et somme le reste
        top_19 = sorted_count[:19]
        others = sorted_count[19:]
        sum_others = sum(count for name, count in others)

        # 3. Écriture des résultats avec pourcentages
        with open(EXIT_FOLDER, 'w', encoding='utf-8') as f_out:
            f_out.write("Phylum;Reads;Abundance_Percent\n")
            # Écriture du Top 19
            for name, count in top_19:
                percent = (count / total_reads) * 100
                f_out.write(f"{name};{count};{percent:.2f}%\n")
            
            # Écriture de la ligne "Others" (si elle n'est pas vide)
            if sum_others > 0:
                percent_others = (sum_others / total_reads) * 100
                f_out.write(f"Others;{sum_others};{percent_others:.2f}%\n")
        print(f"✅ Terminé : {total_reads} reads analysés.")
        return True, total_reads
    return False

# --- EXÉCUTION (À la fin du script) ---
if __name__ == "__main__":
    # Supprime physiquement le fichier avant de le recréer
    if os.path.exists(EXIT_FOLDER):
        shutil.rmtree(EXIT_FOLDER)
        os.makedirs(EXIT_FOLDER)
        print(f"Mise à jour du dossier {EXIT_FOLDER}.")

    # Parcours des files
    for files in os.listdir(ROOT_FOLDER):
        if files.endswith(".kaijuNR"):
            path_in = os.path.join(ROOT_FOLDER, files)
            name_exit = f"abundance_phylum_{files.replace('.kaijuNR', '.txt')}"
            path_out = os.path.join(EXIT_FOLDER, name_exit)
            
            print(f"Analyse de {files} en cours...")
            succes = phylum_abundance(ZIP_PATH, NAME_FILE_IN_ZIP, path_in, path_out)
            
            if succes:
                print(f"✅ Terminé : {name_exit}")
            else:
                print(f"⚠️ Aucun read trouvé dans {files}")
print(f"✅ Analyse terminée")