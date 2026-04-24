import pandas as pd
import os
import sys

def annotate_with_hierarchy(input_dir, ref_path, output_dir):
    """
    Maps KO counts to KEGG hierarchy and applies weight-based abundance adjustment.
    Optimized for memory efficiency and robust metadata handling.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 1. Load and Optimize Reference
    try:
        ref_df = pd.read_csv(ref_path, sep='\t')
        # Clean column names and ensure KO format is consistent
        ref_df.columns = ref_df.columns.str.strip()
        ref_df['ko'] = ref_df['ko'].astype(str).str.replace('ko:', '', regex=False)
        
        # Optimization: Pre-fill NaN values in the reference to avoid doing it in the loop
        ref_df[['level_1', 'level_2', 'level_3']] = ref_df[['level_1', 'level_2', 'level_3']].fillna('Unassigned')
        ref_df['weight'] = ref_df['weight'].fillna(0)
    except Exception as e:
        print(f"Error loading reference file: {e}", file=sys.stderr)
        return

    # 2. Process Samples
    files = [f for f in os.listdir(input_dir) if f.endswith('.tsv')]
    if not files:
        print(f"No files found in {input_dir}")
        return

    for filename in files:
        path = os.path.join(input_dir, filename)
        try:
            # Load sample data
            sample_df = pd.read_csv(path, sep='\t')
            
            # Standardization of join column
            if 'kegg' in sample_df.columns:
                sample_df.rename(columns={'kegg': 'ko'}, inplace=True)
            
            sample_df['ko'] = sample_df['ko'].astype(str).str.replace('ko:', '', regex=False)

            # --- OPTIMIZED MERGE ---
            # Perform a left join to keep all sample KOs, even those without hierarchy info
            merged_df = sample_df.merge(ref_df, on="ko", how="left")

            # Identify numeric columns for adjustment
            # We exclude 'weight' and 'ko' (if interpreted as numeric)
            cols_to_adjust = sample_df.select_dtypes(include=['number']).columns
            
            # --- VECTORIZED WEIGHT ADJUSTMENT ---
            for col in cols_to_adjust:
                # Apply weight adjustment: Adjusted_Count = Original_Count * Weight
                # This correctly splits abundance across multiple pathways
                merged_df[f'adj_{col}'] = merged_df[col] * merged_df['weight']
                
                # Option: Drop the original column to save space if needed
                # merged_df.drop(columns=[col], inplace=True)

            # Final cleanup of missing values for unannotated KOs
            fill_values = {
                'level_1': 'Unassigned',
                'level_2': 'Unassigned',
                'level_3': 'Unassigned',
                'weight': 0
            }
            merged_df.fillna(value=fill_values, inplace=True)

            # Save result
            output_path = os.path.join(output_dir, f"annotated_{filename}")
            merged_df.to_csv(output_path, sep="\t", index=False)
            
            print(f"Successfully annotated: {filename}")

        except Exception as e:
            print(f"Error processing {filename}: {e}", file=sys.stderr)

# --- Execution ---
if __name__ == "__main__":
    annotate_with_hierarchy(
        input_dir='Statistics/35_KO_agreg/',
        ref_path='Statistics/pathway/pathway_levels_extracted.tsv',
        output_dir='Statistics/36_Kegg_merge_pathway_levels/'
    )