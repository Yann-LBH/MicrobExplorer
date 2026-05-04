import csv
from types import SimpleNamespace

# 1. --- COMPATIBILITY SETTINGS ---
try:
    snakemake
except NameError:
    snakemake = SimpleNamespace(
        input=SimpleNamespace(), 
        output=SimpleNamespace()
    )

# 2. --- JOB DESCRIPTION ---
def get_union_and_extract(path_abund, path_rpkm_filt, path_source, path_out):
    """
    1. Collects unique IDs from both filtered files (Union).
    2. Extracts full data for these IDs from the source RPKM file.
    """
    target_ids = set()

    try:  
        # Step 1: Create the Union of IDs
        # Reading filtered abundance
        with open(path_abund, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for cols in reader:
                if cols and not cols[0].lower().startswith("Contig_ID"):
                    target_ids.add(cols[0])
        
        # Adding IDs from filtered RPKM
        with open(path_rpkm_filt, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            for cols in reader:
                if cols and not cols[0].lower().startswith("Contig_ID"):
                    target_ids.add(cols[0])

        # Step 2: Extract data from the source file
        count_final = 0
        with open(path_source, 'r', encoding='utf-8') as f_in, \
            open(path_out, 'w', encoding='utf-8', newline='') as f_out:
            
            reader = csv.reader(f_in, delimiter='\t')
            writer = csv.writer(f_out, delimiter='\t')
            
            for cols in reader:
                # Always keep and write the header
                if cols and cols[0].lower().startswith("Contig_ID"):
                    writer.writerow(cols)
                    continue
                
                # Write row if ID is in our union set
                if cols and cols[0] in target_ids:
                    writer.writerow(cols)
                    count_final += 1
                        
        return count_final

    except Exception as e:
        print(f"❌ Error during extraction: {e}")
        return None

# --- VARIABLE RETRIEVAL & EXECUTION ---
if __name__ == "__main__":

path_abund = snakemake.input.raw_data[0]
current_file = [f for f in path_abund if snakemake.wildcards.sample in f][0]
path_rpkm_filt = snakemake.input.raw_data[1]  
path_out = snakemake.output.rpkm

print(f"🔄 Processing union for output: {snakemake.output.rpkm}")
# Running the intersection logic
if get_union_and_extract(current_file, path_abund, path_rpkm_filt, path_out):
    print(f"✅ Intersection successful -> {path_out}")
else:
    # Raising an error ensures Snakemake identifies the job as failed
    raise RuntimeError(f"Intersection failed for {current_file}")
