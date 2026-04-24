import re
import sys
from collections import defaultdict

def verify_kegg_duplicates(gff_file, feature_type='CDS'):
    """
    Checks for multiple distinct KEGG annotations on the same genomic region.
    Optimized with sets for O(1) lookups and memory efficiency.
    """
    # Dictionary structure: {(contig, start, end): {set_of_kegg_ids}}
    regions = defaultdict(set)
    kegg_pattern = re.compile(r'K\d{5}')
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                # Fast skip for metadata
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.rstrip().split('\t')
                if len(parts) < 9:
                    continue
                
                # Filter by feature type (e.g., CDS) early to save processing
                if parts[2] != feature_type:
                    continue
                
                contig, start, end, attributes = parts[0], parts[3], parts[4], parts[8]
                
                kegg_matches = kegg_pattern.findall(attributes)
                
                if kegg_matches:
                    region_key = (contig, start, end)
                    # Using set.update() is faster than a manual loop with 'if not in'
                    regions[region_key].update(kegg_matches)

        # Analysis of results
        duplicate_count = 0
        print(f"--- KEGG Annotation QC: {gff_file} ---")
        
        for (contig, start, end), keggs in regions.items():
            if len(keggs) > 1:
                duplicate_count += 1
                # Use a single print with join for efficiency
                kegg_list = ", ".join(sorted(keggs))
                print(f"ALERT: Region {contig} [{start}-{end}] has multiple KEGGs: {kegg_list}")
                
        if duplicate_count == 0:
            print("Success: No regions found with multiple distinct KEGG IDs.")
        else:
            print(f"Total: {duplicate_count} regions with multiple annotations found.")

    except FileNotFoundError:
        print(f"Error: File {gff_file} not found.", file=sys.stderr)

# --- Execution ---
if __name__ == "__main__":
    # You can easily change the file or feature type here
    verify_kegg_duplicates('data/coassembly_result.gff3', feature_type='CDS')