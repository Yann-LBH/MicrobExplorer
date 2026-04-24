import pandas as pd
import glob
import os
import sys

def intersection_kegg(fichier_kegg_tsv, dossier_counts, dossier_sortie, extension=".tsv"):
    """
    Optimized intersection between KEGG annotations and count files.
    Uses indexing for faster joins and handles large datasets efficiently.
    """
    if not os.path.exists(dossier_sortie):
        os.makedirs(dossier_sortie)

    # 1. Load and Optimize KEGG Reference
    try:
        # Load only necessary columns if your KEGG file is very large
        df_kegg = pd.read_csv(fichier_kegg_tsv, sep='\t')
        
        # Standardize join column: strip spaces and force string type
        df_kegg['contig'] = df_kegg['contig'].astype(str).str.strip()
        
        # Setting 'contig' as index significantly speeds up the merge operation
        df_kegg.set_index('contig', inplace=True)
    except Exception as e:
        print(f"Error loading KEGG file: {e}")
        return

    # 2. Identify target files
    fichiers_counts = glob.glob(os.path.join(dossier_counts, f"*{extension}"))
    if not fichiers_counts:
        print(f"No files found in {dossier_counts} with extension {extension}")
        return

    print(f"Processing {len(fichiers_counts)} files...")

    # 3. Processing Loop
    for f in fichiers_counts:
        nom_fichier = os.path.basename(f)
        
        try:
            # Load count data
            df_temp = pd.read_csv(f, sep='\t')
            
            # Ensure the join column exists and is clean
            if 'Contig_ID' not in df_temp.columns:
                print(f"Warning: 'Contig_ID' not found in {nom_fichier}. Skipping.")
                continue
                
            df_temp['Contig_ID'] = df_temp['Contig_ID'].astype(str).str.strip()

            # Clean unnecessary columns using a set intersection for speed
            cols_to_drop = {'Reads_Unmapped'}
            existing_drops = cols_to_drop.intersection(df_temp.columns)
            if existing_drops:
                df_temp.drop(columns=list(existing_drops), inplace=True)

            # 4. Optimized Join
            # Using 'left_on' with a pre-indexed 'right_index' is extremely fast
            df_intersection = df_temp.merge(
                df_kegg, 
                left_on='Contig_ID', 
                right_index=True, 
                how='inner'
            )

            if df_intersection.empty:
                print(f"Note: No common contigs found for {nom_fichier}")
                continue

            # 5. Save Output
            clean_name = nom_fichier.replace('Raw_counting_', '')
            chemin_sortie = os.path.join(dossier_sortie, f"kegg_intersect_{clean_name}")
            
            # Export with compression if disk space is an issue (optional)
            df_intersection.to_csv(chemin_sortie, sep='\t', index=False)
            print(f"Successfully generated: {chemin_sortie}")

        except Exception as e:
            print(f"Failed to process {nom_fichier}: {e}")

# --- Execution ---
if __name__ == "__main__":
    intersection_kegg(
        'Statistics/1_Kegg_filtered/kegg_filtered.tsv', 
        'data/count_contigs/', 
        'Statistics/2_Intersection/'
    )