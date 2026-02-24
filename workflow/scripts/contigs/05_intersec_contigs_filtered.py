import os

# --- CONFIGURATION ---
FOLDER_FILTERED = "Statistics/2.Filtered_Contigs"
FOLDER_RPKM = "Statistics/4.RPKM_Contigs_Trimmed"
EXIT_FOLDER = "Statistics/5.Final_Intersection"

def intersect_filters(path_abund, path_rpkm, path_out):
    """ Garde les lignes du fichier RPKM seulement si le contig est présent dans le fichier Filtered_Contigs. """
    try:
        # 1. Charger les IDs valides du fichier Filtered_Contigs dans un set
        valid_ids = set()
        with open(path_abund, 'r', encoding='utf-8') as f:
            for line in f:
                cols = line.strip().split('\t')
                if cols:
                    valid_ids.add(cols[0])

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
        return count_final
    except Exception as e:
        print(f"❌ Erreur sur l'intersection : {e}")
        return None

if __name__ == "__main__":
    os.makedirs(EXIT_FOLDER, exist_ok=True)

    # On liste les fichiers du dossier RPKM pour servir de base
    for filename in os.listdir(FOLDER_RPKM):
        if filename.endswith(".txt"):
            # On cherche le fichier correspondant dans le dossier Abondance
            # On adapte le nom si besoin (ex: RPKM_TD1.txt vs Abundance_filtered_TD1.txt)
            sample_id = filename.replace("RPKM_Contigs_Trimmed_", "").replace(".txt", "")
            
            # Construction des chemins (attention aux noms exacts de tes fichiers)
            path_rpkm = os.path.join(FOLDER_RPKM, filename)
            # Ici, on suppose que le fichier abondance s'appelle ainsi :
            path_abund = os.path.join(FOLDER_FILTERED, f"Filtered_Contigs_{sample_id}.txt")
            
            path_out = os.path.join(EXIT_FOLDER, f"Final_Selected_{sample_id}.txt")

            if os.path.exists(path_abund):
                print(f"🔄 Intersection pour l'échantillon {sample_id}...")
                result = intersect_filters(path_abund, path_rpkm, path_out)
                if result is not None:
                    print(f"✅ Terminé : {result} contigs conservés.")
            else:
                print(f"⚠️ Fichier filtered introuvable pour {sample_id}")

    print(f"\n🚀 Dossier final créé : {EXIT_FOLDER}")