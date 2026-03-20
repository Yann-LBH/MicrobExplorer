import csv
from types import SimpleNamespace

# 1. --- COMPATIBILITY SETTINGS ---
try:
    snakemake
except NameError:
    snakemake = SimpleNamespace(
        input=SimpleNamespace(), 
        output=SimpleNamespace(),
        wildcards=SimpleNamespace()
    )

# 2. --- JOB DESCRIPTION ---
def load_taxonomy(tax_path):
    """
    Loads taxonomy TSV into a dictionary {Contig_ID: Lineage}.
    Expected format: [0] Class | [1] ID | [2] TaxID | [3] Lineage
    """
    tax_dict = {}
    try:
        with open(tax_path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for cols in reader:
                if len(cols) >= 4:
                    contig_id = cols[1].strip()
                    lineage = cols[3].strip()
                    # Clean trailing semicolon if present
                    if lineage.endswith(';'):
                        lineage = lineage[:-1]
                    tax_dict[contig_id] = lineage
    except FileNotFoundError:
        print(f"⚠️ Warning: Taxonomy file not found at {tax_path}")
    return tax_dict

def run_annotation(input_path, output_path, tax_dict):
    """
    Annotates RPKM data with taxonomy and sorts results by RPKM descending.
    """
    data_to_sort = []
    header = None

    try:
        with open(input_path, 'r', encoding='utf-8') as f_in:
            reader = csv.reader(f_in, delimiter='\t')
            
            # Capture header
            try:
                header = next(reader)
                header.append("Taxonomy")
            except StopIteration:
                return False

            for cols in reader:
                if len(cols) < 4:
                    continue
                
                contig_id = cols[0].strip()
                lineage = tax_dict.get(contig_id, "Unclassified")
                
                try:
                    # Column index 3 is RPKM
                    rpkm_val = float(cols[3])
                    # Add lineage to the row data
                    cols.append(lineage)
                    data_to_sort.append((rpkm_val, cols))
                    
                except ValueError:
                    continue

        # Sort data by RPKM (index 0 of the tuple) in descending order
        data_to_sort.sort(key=lambda x: x[0], reverse=True)

        # Write annotated and sorted file
        with open(output_path, 'w', encoding='utf-8', newline='') as f_out:
            writer = csv.writer(f_out, delimiter='\t')
            if header:
                writer.writerow(header)
            
            for item in data_to_sort:
                writer.writerow(item[1])

        return True

    except Exception as e:
        print(f"❌ Error during annotation: {e}")
        return False

# --- VARIABLE RETRIEVAL & EXECUTION ---

path_in = snakemake.input.data
current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
path_taxonomy = snakemake.input.taxonomy
path_out = snakemake.output.taxaname

print(f"📊 Loading taxonomy reference...")
taxonomy_map = load_taxonomy(path_taxonomy)

print(f"📊 Annotating {current_file}...")
if run_annotation(current_file, path_out, taxonomy_map):
    print(f"✅ Annotation successful -> {path_out}")
else:
    raise RuntimeError(f"Annotation failed for {current_file}")