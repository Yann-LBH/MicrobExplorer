import os
import re

# --- CONFIGURATION DES CHEMINS ---
STEPS = {
    "0_Brut": "Data",               # .tsv
    "1_Counted": "Statistics/1.Raw_Counting",  # .txt
    "2_Trimmed": "Statistics/2.Filtered_Contigs", 
    "3_RPKM": "Statistics/3.RPKM_Contigs",
    "4_RPKM_Filtered": "Statistics/4.RPKM_Contigs_Trimmed",
    "5_Intersection": "Statistics/5.Intersection", # .csv ou .txt
    "6_Final": "Statistics/6.Final_Results_CSV" # .csv ou .txt
}

def get_contig_ids(path, is_csv=False, col_index=0):
    """Extrait l'ensemble des IDs uniques de contigs d'un fichier."""
    ids = set()
    if not os.path.exists(path): return ids
    
    with open(path, 'r', encoding='utf-8') as f:
        # Si c'est un TSV/CSV, on ignore l'en-tête
        if is_csv: next(f, None) 
        
        for line in f:
            line = line.strip()
            if not line: continue
            
            # Gestion du séparateur : tab pour tsv, ; pour csv, sinon auto
            if path.endswith('.tsv'):
                parts = line.split('\t')
            elif is_csv:
                parts = line.split(';')
            else:
                parts = line.split() # Par défaut : espaces/tabs
            
            if parts and len(parts) > col_index:
                # Nettoyage des guillemets et espaces
                contig_id = parts[col_index].strip('"').strip()
                if contig_id: ids.add(contig_id)
    return ids

def run_full_qc():
    qc_results = {}
    all_samples = set()

    # --- ÉTAPE 1 À 6 : ANALYSE ---
    for step_key, folder in STEPS.items():
        print(f"Analyse {step_key}...")
        
        if not os.path.exists(folder): 
            print(f"  ⚠️ Dossier {folder} introuvable. Skip.")
            continue
            
        for f in os.listdir(folder):
            # 1. Nettoyage du nom de l'échantillon
            name = f.replace(".tsv", "").replace(".csv", "")
            prefixes = ["count-contigs-coassembly-", "Raw_counting_", "Filtered_Contigs_", "RPKM_Contigs_", "Trimmed_", "Final_Selected_"]
            #print(f"DEBUG: Fichier '{f}' classé comme Sample '{name}'")
            for prefix in prefixes:
                name = name.replace(prefix, "")
            
            all_samples.add(name)
            
            # 2. Détection du format
            is_csv_or_tsv = f.endswith(".csv") or f.endswith(".tsv")
            path = os.path.join(folder, f)
            
            # 3. Récupération des IDs
            ids = get_contig_ids(path, is_csv=is_csv_or_tsv)
            
            if name not in qc_results:
                qc_results[name] = {}
            
            # On stocke le nombre pour le rapport
            qc_results[name][step_key] = len(ids)

    # --- AFFICHAGE DU RAPPORT ---
    print("\n" + "="*100)
    header = f"{'Sample':<20} | {'Brut':<8} | {'Counted':<8} | {'Trimmed':<8} | {'RPKM':<8} | {'RPKM_Filtered':<8} | {'Intersec':<8} | Perte Total | Statut"
    print(header)
    print("-" * 100)

    for sample in sorted(all_samples):
        res = qc_results.get(sample, {})
        
        # Ici, on récupère le NOMBRE de contigs restants à chaque dossier
        v0 = res.get("0_Brut", 0)           # Stock initial
        v1 = res.get("1_Counted", 0)        # Après premier comptage
        v2 = res.get("2_Trimmed", 0)        # Après filtre de longueur
        v3 = res.get("3_RPKM", 0)           # Après calcul RPKM
        v4 = res.get("4_RPKM_Filtered", 0)  # Après filtre RPKM > 1
        v5 = res.get("5_Intersection", 0)   # Stock final après intersection
        v6 = res.get("6_Final", 0)
        
        pertes_totales = (v0 - v1) + (v1 - v2) + (v2 - v3) + (v3 - v4) + (v4 - v5)
        status = "✅ OK" if pertes_totales + v5 == v0 else "⚠️ ERROR"

        print(f"{sample:<20} | {v1:<8} | {v2:<8} | {v3:<8} | {v4:<8} | {v5:<8} | {v6:<8} | {pertes_totales:<8} | {status}")

    print("="*100)

if __name__ == "__main__":
    run_full_qc()