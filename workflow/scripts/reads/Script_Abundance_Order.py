import os
import zipfile
import shutil

# --- CONFIGURATION (En haut pour être facilement modifiable) ---
ZIP_PATH = "A.Taxa_Table_NCBI/new_taxdump.zip"
NAME_FILE_IN_ZIP = "rankedlineage.dmp"
ROOT_FOLDER = "Data"
EXIT_FOLDER = "Abundance/Order"

def order_abundance(ZIP_PATH, NAME_FILE_IN_ZIP, ROOT_FOLDER, EXIT_FOLDER):
    # 1. Create output directory if it doesn't exist    
    output_dir = os.path.dirname(EXIT_FOLDER)

    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"📁 Dossier créé : {output_dir}")
        
    print(f"Détermination des ordres des reads depuis {NAME_FILE_IN_ZIP}...")
    tax_to_order = {}
    
    with zipfile.ZipFile(ZIP_PATH, 'r') as z:
        with z.open(NAME_FILE_IN_ZIP) as f_in:
            
            for binary_line in f_in:
                parts = [p.strip() for p in binary_line.decode('utf-8').split('|')]

                if len(parts) >= 9:
                    tax_to_order[parts[0]] = parts[5] if parts[5] else "Unknown Order" #Associe TaxaId au Domain 0 = tax_id, 9 = domain
    
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
                order = tax_to_order.get(tax_id, "Unknown/Other") #Si le tax_id n'est pas dans le mapping
                counter[order] = counter.get(order, 0) + 1
            else:
                counter["Unclassified"] = counter.get("Unclassified", 0) + 1

    # 3. Écriture des résultats avec pourcentages
    if total_reads > 0:
        with open(EXIT_FOLDER, 'w', encoding='utf-8') as f_out:
            f_out.write("Order;Reads;Abundance_Percent\n")
            for order, count in sorted(counter.items(), key=lambda x: x[1], reverse=True):
                percent = (count / total_reads) * 100
                f_out.write(f"{order};{count};{percent:.2f}%\n")
        print(f"✅ Terminé : {total_reads} reads analysés.")
        return True
    return False

# --- EXÉCUTION (À la fin du script) ---
if __name__ == "__main__":
    # Supprime physiquement le fichier avant de le recréer
    if os.path.exists(EXIT_FOLDER):
        shutil.rmtree(EXIT_FOLDER)
        os.makedirs(EXIT_FOLDER)
        print(f"Mise à jour du dossier {EXIT_FOLDER}.")

    # Parcours des filess
    for files in os.listdir(ROOT_FOLDER):
        if files.endswith(".kaijuNR"):
            path_in = os.path.join(ROOT_FOLDER, files)
            name_exit = f"abundance_order_{files.replace('.kaijuNR', '.txt')}"
            path_out = os.path.join(EXIT_FOLDER, name_exit)
            
            print(f"Analyse de {files} en cours...")
            succes = order_abundance(ZIP_PATH, NAME_FILE_IN_ZIP, path_in, path_out)
            
            if succes:
                print(f"✅ Terminé : {name_exit}")
            else:
                print(f"⚠️ Aucun read trouvé dans {files}")
print(f"✅ Analyse terminée")