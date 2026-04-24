import re
import csv
import sys
import argparse

def process_gff_kegg(input_file, output_file):
    """
    Parses a GFF3 file and extracts KEGG annotations using streaming to save memory.
    """
    # Pre-compile regex for speed
    kegg_pattern = re.compile(r'K\d{5}')
    
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w', newline='') as f_out:
            writer = None
            count = 0
            
            for line in f_in:
                # Fast skip for comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                
                parts = line.rstrip().split('\t')
                if len(parts) < 9:
                    continue
                
                # Column assignment
                contig = parts[0]
                try:
                    # GFF coordinates are 1-based inclusive
                    start = int(parts[3])
                    end = int(parts[4])
                    length = end - start + 1
                except ValueError:
                    continue

                attributes = parts[8]
                keggs_found = set(kegg_pattern.findall(attributes))
                
                for kegg in keggs_found:
                    row = {
                        'contig': contig,
                        'kegg': kegg,
                        'gene_length': length
                    }
                    
                    # Initialize CSV writer with the first valid row found
                    if writer is None:
                        writer = csv.DictWriter(f_out, fieldnames=row.keys(), delimiter='\t')
                        writer.writeheader()
                    
                    writer.writerow(row)
                    count += 1
            
            print(f"Extraction completed: {count} annotations processed.")

    except FileNotFoundError:
        print(f"Error: File {input_file} not found.", file=sys.stderr)
    except PermissionError:
        print(f"Error: Permission denied for output file.", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Extract KEGG annotations from GFF3 files.")
    parser.add_argument("-i", "--input", required=True, help="Path to input GFF3 file")
    parser.add_argument("-o", "--output", required=True, help="Path to output TSV file")
    
    args = parser.parse_args()
    process_gff_kegg(args.input, args.output)

if __name__ == "__main__":
    main()