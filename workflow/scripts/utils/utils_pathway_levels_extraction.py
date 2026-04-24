import requests
import pandas as pd
import os
import re
import sys

def get_kegg_hierarchy():
    """
    Fetches and parses the KEGG hierarchy (ko00001) to extract 
    functional levels, EC numbers, and calculates KO weights.
    """
    url = "https://rest.kegg.jp/get/br:ko00001"
    
    # Pre-compile regex patterns for significant speedup in loops
    re_clean_ab = re.compile(r'^[A-B\s]*\d+\s+')
    re_pathway = re.compile(r'^\d+\s+(.*)\s*\[PATH:ko\d+\]')
    re_ec = re.compile(r'\[EC:(.*?)\]')
    
    try:
        # Increased timeout and robust request handling
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        hierarchy_data = []
        l1, l2, l3 = "", "", ""
        
        # Split text into lines once to avoid repeated split calls
        lines = response.text.splitlines()

        for line in lines:
            if not line:
                continue
                
            prefix = line[0]
            
            if prefix == 'A':
                l1 = re_clean_ab.sub('', line).strip()
            elif prefix == 'B':
                l2 = re_clean_ab.sub('', line).strip()
            elif prefix == 'C':
                # Fast extraction of pathway name
                match_c = re_pathway.search(line[3:])
                l3 = match_c.group(1).strip() if match_c else line[3:].strip()
            elif prefix == 'D':
                content = line[4:].strip()
                
                # Extract EC numbers if they exist
                ec_match = re_ec.search(content)
                ec_val = ec_match.group(1) if ec_match else ""
                
                # Clean content for ID and Description
                clean_content = re_ec.sub('', content).strip()
                
                # Use maxsplit=1 for efficiency if we only care about the first ID
                parts = clean_content.split(None, 1)
                ko_id = parts[0]
                ko_desc = parts[1] if len(parts) > 1 else ""
                
                hierarchy_data.append({
                    "ko": ko_id,
                    "ec_number": ec_val,
                    "level_1": l1,
                    "level_2": l2,
                    "level_3": l3,
                    "gene_description": ko_desc
                })

        # --- DataFrame Optimization ---
        df = pd.DataFrame(hierarchy_data)
        
        # Calculate weight: 1 / frequency of KO across all pathways
        # This prevents overestimating the abundance of multi-functional KOs
        df['weight'] = 1 / df.groupby('ko')['ko'].transform('count')
        
        return df

    except requests.exceptions.RequestException as e:
        print(f"Network error fetching KEGG data: {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Parsing error: {e}", file=sys.stderr)
        return None

def main():
    output_dir = "Statistics/pathway"
    output_file = "pathway_levels_extracted.tsv"
    
    print("Connecting to KEGG database...")
    kegg_df = get_kegg_hierarchy()
    
    if kegg_df is not None:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, output_file)
        
        # Save to TSV
        kegg_df.to_csv(path, sep='\t', index=False)
        
        print(f"Successfully generated hierarchy: {path}")
        print(f"Total entries: {len(kegg_df)}")
        print(kegg_df.head(3))

if __name__ == "__main__":
    main()