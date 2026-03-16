import csv

# --- SNAKEMAKE CONFIGURATION (VSCode alerts) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(raw_data=["filtered_abund.tsv", "sample_rpkm.tsv"]), 
        output=SimpleNamespace(rpkm="filtered_output.tsv")
    )

# --- VARIABLE RETRIEVAL ---
# Accessing the two specific input files from the Snakemake list
path_abund = snakemake.input.raw_data[0]
current_file = [f for f in path_abund if snakemake.wildcards.sample in f][0]
path_rpkm = snakemake.input.raw_data[1]  
path_out = snakemake.output.rpkm

def intersect_filter(path_abund, path_rpkm, path_out):
    """
    Filters the RPKM file by keeping only contigs present in the abundance reference file.
    """
    try:
        # 1. Load reference IDs into a set for O(1) fast lookup
        valid_ids = set()
        with open(path_abund, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)  # Skip reference header
            for cols in reader:
                if cols:
                    # Collect the ID from the first column
                    valid_ids.add(cols[0])

        # 2. Filter the RPKM file based on the collected set
        with open(path_rpkm, 'r', encoding='utf-8') as f_in, \
             open(path_out, 'w', encoding='utf-8', newline='') as f_out:
            
            reader = csv.reader(f_in, delimiter='\t')
            writer = csv.writer(f_out, delimiter='\t')

            # Handle the Header for the output file
            header = next(reader, None)
            if header:
                writer.writerow(header)

            # Filtering rows and counting preserved entries
            count_final = 0
            for cols in reader:
                if cols and cols[0] in valid_ids:
                    writer.writerow(cols)
                    count_final += 1
            
            print(f"📈 {count_final} contigs preserved in the final intersection.")
        
        return True

    except Exception as e:
        print(f"❌ Error during intersection process: {e}")
        return False

# --- EXECUTION ---
# Running the intersection logic
if intersect_filter(current_file, path_abund, path_rpkm, path_out):
    print(f"✅ Intersection successful -> {path_out}")
else:
    # Raising an error ensures Snakemake identifies the job as failed
    raise RuntimeError(f"Intersection failed for {current_file}")