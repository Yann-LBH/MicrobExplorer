import pandas as pd
import glob
import os
import sys

def sum_per_kegg(input_dir, output_dir, extension=".tsv"):
    """
    Aggregates read counts by KEGG ID and prepares the matrix for differential expression.
    Optimized for memory by selecting columns during load and using efficient grouping.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Identify files from the previous intersection step
    files = glob.glob(os.path.join(input_dir, f"*{extension}"))
    
    if not files:
        print(f"No files found in {input_dir}")
        return

    # Metadata columns to exclude from summation
    metadata_cols = {'contig', 'Length', 'Contig_ID', 'gene_length', 'id', 'product'}

    for f in files:
        file_name = os.path.basename(f)
        
        try:
            # 1. Optimization: Read only the header first to filter columns
            # This avoids loading massive metadata columns into RAM
            header = pd.read_csv(f, sep='\t', nrows=0).columns
            cols_to_keep = [c for c in header if c == 'kegg' or c not in metadata_cols]
            
            # 2. Load only the necessary data
            df = pd.read_csv(f, sep='\t', usecols=cols_to_keep)

            if 'kegg' not in df.columns:
                print(f"Warning: 'kegg' column missing in {file_name}. Skipping.")
                continue

            # 3. Aggregation
            # numeric_only=True ensures we don't try to sum strings by mistake
            df_final = df.groupby('kegg').sum(numeric_only=True).reset_index()

            # 4. Save results
            output_name = f"final_counts_{file_name}"
            output_path = os.path.join(output_dir, output_name)
            df_final.to_csv(output_path, sep='\t', index=False)
            
            print(f"Matrix ready for downstream analysis: {output_path}")

        except Exception as e:
            print(f"Error processing {file_name}: {e}", file=sys.stderr)

# --- Execution ---
if __name__ == "__main__":
    sum_per_kegg(
        'Statistics/2_Intersection/', 
        'Statistics/3_Deseq2/'
    )