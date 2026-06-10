################################################################################
# Project : "MicrobExplorer"
# Script: " utils Pathway levels extraction"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

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
    # Retrieve the URL provided by the Snakemake execution context
    url = snakemake.params.url
    
    # Pre-compile regex patterns for significant speedup in loops
    re_clean_ab = re.compile(r'^[A-B\s]*\d+\s+')
    re_pathway  = re.compile(r'^\d+\s+(.*)\s*\[PATH:ko\d+\]')
    re_ec       = re.compile(r'\[EC:(.*?)\]')
    
    # ✅ FIX: All networking logic is now properly wrapped inside the try block
    try:
        print(f"Fetching data from KEGG API: {url}")
        # Single robust request with mandatory timeout
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        hierarchy_data = []
        l1, l2, l3 = "", "", ""
        
        # Split text into lines once to optimize iteration performance
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
                # ✅ FIX: Isolated local variable name to prevent confusion
                line_content = line[4:].strip()
                
                # Extract EC numbers if they exist
                ec_match = re_ec.search(line_content)
                ec_val = ec_match.group(1) if ec_match else ""
                
                # Clean content for ID and Description extraction
                clean_content = re_ec.sub('', line_content).strip()
                
                # Efficient split to capture KO identifier and its associated description
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

        # --- DataFrame Optimization & Scaling ---
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
    # Dynamic path fetching from snakemake output if available, otherwise defaults to local
    output_file = getattr(snakemake.output, "tsv", "Statistics/pathway/pathway_levels_extracted.tsv")
    output_dir = os.path.dirname(output_file)
    
    print("Connecting to KEGG database...")
    kegg_df = get_kegg_hierarchy()
    
    if kegg_df is not None:
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        # Save results directly to the targeted file path
        kegg_df.to_csv(output_file, sep='\t', index=False)
        
        print(f"✓ Successfully generated hierarchy: {output_file}")
        print(f"Total entries: {len(kegg_df)}")
        print(kegg_df.head(3))
    else:
        # Fail explicitly so Snakemake stops the execution flow
        raise RuntimeError("KEGG hierarchy extraction failed. Review error logs above.")

if __name__ == "__main__":
    main()