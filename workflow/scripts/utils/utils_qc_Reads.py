################################################################################
# Project : "MicrobExplorer"
# Script: "Utils for Reads QC"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################
import os
import re

def get_total_reads(path, is_csv=False):
    """Calcule le total des reads d'un fichier en un seul passage."""
    total = 0
    if not os.path.exists(path): 
        return 0
    
    with open(path, 'r', encoding='utf-8') as f:
        if is_csv:
            next(f, None)  # Skip TSV/CSV header
        for line in f:
            if is_csv:
                # Adapte le séparateur si ton fichier final utilise des tabulations ou des points-virgules
                parts = line.strip().split('\t' if '\t' in line else ';')
                if len(parts) >= 2:
                    val = re.sub(r'[^0-9]', '', parts[1])
                    if val: 
                        total += int(val)
            else:
                # Find the last sequence of digits in the line (reads count)
                numbers = re.findall(r'\d+', line)
                if numbers:
                    total += int(numbers[-1])
    return total

def run_full_qc(snakemake):
    """
    Main execution function for Reads QC.
    Processes the exact files provided by the Snakemake input list.
    """
    qc_results = {}

    # 1. ÉTAPE 0 : COMPTAGE KAIJU BRUT (.kaijuNR)
    print("Analyzing Step 0 (Raw Data)...")
    brut_files = [f for f in DATA if f.endswith(".kaijuNR")]
    for path in brut_files:
        filename = os.path.basename(path)
        sample_name = filename.replace(".kaijuNR", "").replace("reads_", "")
        
        with open(path, 'r', encoding='utf-8') as file:
            count = sum(1 for _ in file) # underscore for unused variable, just counting lines
        qc_results[sample_name] = {"Brut": count}

    # 2. ÉTAPE 1, 2 & 3 : COMPTAGE TEXTE ET FINAL (.tsv)
    # On isole les fichiers TSV pour les traiter selon leur pattern de nommage
    tsv_files = [f for f in DATA if f.endswith(".tsv")]
    
    for path in tsv_files:
        filename = os.path.basename(path)
        
        # Identification de la catégorie de l'étape par rapport au préfixe du fichier
        if filename.startswith("counted_"):
            sample_name = filename.replace("_reads.tsv", "").replace("counted_", "")
            qc_results.setdefault(sample_name, {})["Counted"] = get_total_reads(path, is_csv=False)
            
        elif filename.startswith("filtered_"):
            sample_name = filename.replace("_reads.tsv", "").replace("filtered_", "")
            qc_results.setdefault(sample_name, {})["Filtered"] = get_total_reads(path, is_csv=False)
            
        elif filename.startswith("annotated_"):
            sample_name = filename.replace("_reads.tsv", "").replace("annotated_", "")
            # Étape finale considérée comme un format tabulaire avec en-tête (is_csv=True)
            qc_results.setdefault(sample_name, {})["Final"] = get_total_reads(path, is_csv=True)

    return qc_results

if __name__ == "__main__":
    # Inputs
    DATA = snakemake.input.data
        
    # Ouputs
    PDF     = snakemake.output.pdf
    PARQUET = snakemake.output.parquet
    
    results=run_full_qc(snakemake)
    print(f"✓ Reads QC complete. Data to {PARQUET} for {len(results)} samples.")