import os
import re

# --- CONFIGURATION DES CHEMINS ---
STEPS = {
    "1_Brut": "Data",                # .kaijuNR
    "2_Counted": "1.Counted_kaiju",  # .txt (Taxon X : Y reads)
    "3_Trimmed": "2.Reads_trimmed",  # .txt (Taxon X : Y reads >= 100)
    "4_Final": "3.Finals_Results_CSV" # .csv final
}

def get_total_reads(path, is_csv=False):
    """Calcule le total des reads d'un fichier en un seul passage."""
    total = 0
    if not os.path.exists(path): return 0
    
    with open(path, 'r', encoding='utf-8') as f:
        if is_csv:
            next(f, None) # Sauter l'en-tête
        for line in f:
            if is_csv:
                parts = line.strip().split(';')
                if len(parts) >= 2:
                    # On nettoie les guillemets éventuels et on convertit
                    val = re.sub(r'[^0-9]', '', parts[1])
                    if val: total += int(val)
            else:
                # Cherche le dernier nombre de la ligne (reads)
                numbers = re.findall(r'\d+', line)
                if numbers:
                    total += int(numbers[-1])
    return total

def run_full_qc():
    qc_results = {}

    # 1. ÉTAPE 1 : COMPTAGE KAIJU BRUT (Lignes C)
    print("Analyse Etape 1 (Brut)...")
    for f in os.listdir(STEPS["1_Brut"]):
        if f.endswith(".kaijuNR"):
            name = f.replace(".kaijuNR", "").replace("reads_", "")
            path = os.path.join(STEPS["1_Brut"], f)
            with open(path, 'r') as file:
                count = sum(1 for line in file)
            qc_results[name] = {"Brut": count}

    # 2. ÉTAPE 2 & 3 : COMPTAGE TEXTE (Somme des reads)
    for label, folder, prefix in [("Counted", "2_Counted", "counted_reads_"), ("Trimmed", "3_Trimmed", "trimmed_reads_counted_reads_")]:
        print(f"Analyse Etape {label}...")
        if not os.path.exists(STEPS[folder]): continue
        for f in os.listdir(STEPS[folder]):
            if f.endswith(".txt"):
                name = f.replace(".txt", "").replace(prefix, "")
                qc_results.setdefault(name, {})[label] = get_total_reads(os.path.join(STEPS[folder], f), is_csv=False)

    # 3. ÉTAPE 4 : COMPTAGE CSV (Somme des reads)
    print("Analyse Etape 4 (Final)...")
    if os.path.exists(STEPS["4_Final"]):
        for f in os.listdir(STEPS["4_Final"]):
            if f.endswith(".csv"):
                name = f.replace(".csv", "").replace("final_", "")
                qc_results.setdefault(name, {})["Final"] = get_total_reads(os.path.join(STEPS["4_Final"], f), is_csv=True)

    # --- AFFICHAGE DU RAPPORT ---
    print("\n" + "="*90)
    header = f"{'Sample':<20} | {'Brut':<10} | {'Counted':<10} | {'Trimmed':<10} | {'Final':<10} | {'Statut'}"
    print(header)
    print("-" * 90)
            
    for sample, data in qc_results.items():
        # Récupération sécurisée des valeurs (0 par défaut)
        brut = data.get("Brut", 0)
        cnt  = data.get("Counted", 0)
        trm  = data.get("Trimmed", 0)
        fnl  = data.get("Final", 0)
                
        # Calculs des pertes
        diff_counted = brut - cnt
        loose_trim = cnt - trm
        loose_noise = trm - fnl
                    
        # Statut visuel
        status = "✅ OK" if (fnl + loose_noise + loose_trim + diff_counted) == brut else "⚠️ LOOSE"
                
        print(f"{sample:<20} | {brut:<10,d} | {cnt:<10,d} | {trm:<10,d} | {fnl:<10,d} | {status}")
        if status == "⚠️ LOOSE":
            print(f"    └─ Détails : Filter Counted (-{diff_counted:,d}) | Filter Trimmed (-{loose_trim:,d}) | Filter Final (-{loose_noise:,d})\n")

    print("="*90)

    # --- EXPORT DU RAPPORT DANS UN FICHIER ---
    output_file = "Rapport_QC_Final.txt"
    with open(output_file, 'w', encoding='utf-8') as report:
        report.write("RAPPORT DE QUALITÉ CONTRÔLE (QC)\n")
        report.write("="*90 + "\n")
        report.write(f"{'Sample':<20} | {'Brut':<10} | {'Counted':<10} | {'Trimmed':<10} | {'Final':<10} | {'Statut'}\n")
        report.write("-" * 90 + "\n")
        
        for sample, data in qc_results.items():
            brut = data.get("Brut", 0)
            cnt  = data.get("Counted", 0)
            trm  = data.get("Trimmed", 0)
            fnl  = data.get("Final", 0)
            
            diff_counted = brut - cnt
            loose_trim = cnt - trm
            loose_noise = trm - fnl
            status = "OK" if (fnl + loose_noise + loose_trim + diff_counted) == brut else "LOOSE"
            
            # Écriture de la ligne principale
            report.write(f"{sample:<20} | {brut:<10,d} | {cnt:<10,d} | {trm:<10,d} | {fnl:<10,d} | {status}\n")
            
            # Écriture du détail si anomalie ou perte notable
            if brut != diff_counted + cnt or cnt != loose_trim + trm or trm != loose_noise + fnl:
                report.write(f"    └─ Details : Filter Counted (-{diff_counted:,d}) | Filter Trimmed (-{loose_trim:,d}) | Filter Final (-{loose_noise:,d})\n")
                report.write(f"\n")
        report.write("="*90 + "\n")
    
    print(f"\n✅ Rapport sauvegardé avec succès dans : {output_file}")

if __name__ == "__main__":
    run_full_qc()