import os
import shutil

# --- CONFIGURATION (En haut pour être facilement modifiable) ---
ROOT_FOLDER = "1.Counted_kaiju"
EXIT_FOLDER = "2.Reads_trimmed"
    
def trimming_kaiju(enter_path, exit_path):
    
    try:
        with open(enter_path, 'r', encoding='utf-8') as f_in, \
             open(exit_path, 'w', encoding='utf-8') as f_out:

            #garde les reads >= 100
            for line in f_in:
            # On découpe la ligne par espaces
                columns = line.strip().split(' ')
            
                if len(columns) >= 4:
                    try:
                        # Conversion du nombre de reads en entier
                        number_reads = int(columns[3]) 
                    
                        # On ne garde que si c'est >= 100
                        if number_reads >= 100:
                            f_out.write(line)
                    except ValueError:
                        # En cas d'en-tête (texte au lieu de nombre), on ignore l'erreur
                        continue
            #garde les 100 premières lignes seulement ( les plus abondants)
            #for i, lines in enumerate(f_in):
                #if i >= 100:
                    #break
                #f_out.write(lines)

        return True
    except Exception as e:
            print(f"❌ Erreur sur le fichier {enter_path}: {e}")
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
        if files.endswith(".txt"):
            path_in = os.path.join(ROOT_FOLDER, files)
            name_exit = f"trimmed_reads_{files}"
            path_out = os.path.join(EXIT_FOLDER, name_exit)
            
            print(f"Traitement de {files} en cours...")
            succes = trimming_kaiju(path_in, path_out)
            
            if succes:
                print(f"✅ Terminé : {name_exit}")