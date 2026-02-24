import os
import zipfile

# --- CONFIGURATION ---
ZIP_PATH = "A.Taxa_Table_NCBI/new_taxdump.zip"
NAME_FILE_IN_ZIP = "rankedlineage.dmp"
EXIT_FOLDER = "B.Taxa_Table_CSV"
EXIT_FILE = os.path.join(EXIT_FOLDER, "mapping_taxons.csv")

def taxon_phylo_extraction_from_zip(ZIP_PATH, NAME_FILE_IN_ZIP, exit):
    print(f"Création de la table de correspondance depuis {NAME_FILE_IN_ZIP}...")


    with zipfile.ZipFile(ZIP_PATH, 'r') as z:
        with z.open(NAME_FILE_IN_ZIP) as f_in, open(exit, 'w', encoding='utf-8') as f_out:

            f_out.write("tax_id;scientific_name;domain;kingdom;phylum;class;order;family;genus;species\n")
            
            for binary_line in f_in:
                line = binary_line.decode('utf-8')
                parts = [p.strip() for p in line.split('|')]

                if len(parts) >= 9:
                    tax_id = parts[0]
                    sci_name = parts[1].replace('"', '')
                    # On récupère les colonnes directement (l'ordre NCBI est inversé)
                    # species[2]
                    # genus[3]
                    # family[4]
                    # order[5]
                    # classe[6]
                    # phylum[7]
                    # kingdom[8]
                    # domain[9]

                    # On récupère les colonnes taxonomiques (index 2 à 9)
                    ranks = [p.replace('"', '') for p in parts[2:10]]
                    
                    # --- LOGIQUE DE REMPLISSAGE ---
                    # On parcourt les rangs de gauche à droite (Species -> Kingdom)
                    # Si une colonne est vide, on lui donne le sci_name
                    for i in range(len(ranks)):
                        if ranks[i] == "":
                            ranks[i] = sci_name
                   
                    ranks_ordered = ranks[::-1]

                    # On assemble : ID + SciName + les rangs modifiés
                    final_line = [tax_id, sci_name] + ranks_ordered
                    f_out.write(";".join(final_line) + "\n")
    
    print(f"✅ Table de correspondance créée dans : {exit}")

if __name__ == "__main__":
    # Création du dossier de exit s'il n'existe pas
    if not os.path.exists(EXIT_FOLDER):
        os.makedirs(EXIT_FOLDER)
        print(f"Dossier '{EXIT_FOLDER}' créé.")

# Vérification que le ZIP existe avant de commencer
    if os.path.exists(ZIP_PATH):
        taxon_phylo_extraction_from_zip(ZIP_PATH, NAME_FILE_IN_ZIP, EXIT_FILE)
    else:
        print(f"❌ Erreur : Le fichier {ZIP_PATH} est introuvable.")