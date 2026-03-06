import os

# --- CONFIGURATION (En haut pour être facilement modifiable) ---
ROOT_FOLDER = "Data"
EXIT_FOLDER = "1.Counted_kaiju"

# Supprime physiquement le fichier avant de le recréer
if os.path.exists(EXIT_FOLDER):
    os.remove(EXIT_FOLDER)
    print(f"Ancien fichier {EXIT_FOLDER} supprimé.")

def analyser_kaiju(enter_path, exit_path):
    counter = {}
    
    with open(enter_path, 'r', encoding='utf-8') as f:
        for line in f:
            # Kaiju sépare généralement les column par des tabulations
            column = line.strip().split('\t')
            
            # On vérifie si la line est classée ('C') et possède bien un ID taxon
            if len(column) >= 3 and column[0] == 'C':
                taxon_id = column[2] # La 3ème colonne (ex: 35786)
                
                counter[taxon_id] = counter.get(taxon_id, 0) + 1

    # Écriture des résultats
    if counter:
        with open(exit_path, 'w', encoding='utf-8') as f_out:
            # On trie par nombre d'occurrences (du plus grand au plus petit)
            for taxon, total in sorted(counter.items(), key=lambda x: x[1], reverse=True):
                f_out.write(f"Taxon {taxon} : {total} reads\n")
        return True
    return False

# --- EXÉCUTION (À la fin du script) ---
if __name__ == "__main__":
    # Création du dossier de sortie s'il n'existe pas
    if not os.path.exists(EXIT_FOLDER):
        os.makedirs(EXIT_FOLDER)
        print(f"Dossier '{EXIT_FOLDER}' créé.")

    # Parcours des filess
    for files in os.listdir(ROOT_FOLDER):
        if files.endswith(".kaijuNR"):
            path_in = os.path.join(ROOT_FOLDER, files)
            name_exit = f"counted_{files.replace('.kaijuNR', '.txt')}"
            path_out = os.path.join(EXIT_FOLDER, name_exit)
            
            print(f"Analyse de {files} en cours...")
            succes = analyser_kaiju(path_in, path_out)
            
            if succes:
                print(f"✅ Terminé : {name_exit}")
            else:
                print(f"⚠️ Aucun read classé ('C') trouvé dans {files}")