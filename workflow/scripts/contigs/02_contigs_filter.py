import csv
from types import SimpleNamespace

# 1. --- COMPATIBILITY SETTINGS ---
try:
    snakemake
except NameError:
    snakemake = SimpleNamespace(
        input=SimpleNamespace(), 
        output=SimpleNamespace(), 
        params=SimpleNamespace(), 
        wildcards=SimpleNamespace()
    )

# 2. --- JOB DESCRIPTION ---
def get_total_abundance(files):
    """
    Scan all files to sum reads per contig ID.
    """
    global_counts = {}
    for path in files:
        with open(path, 'r', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            next(reader, None)  # Skip header
            for columns in reader:
                if len(columns) >= 3:
                    try:
                        c_id, reads = columns[0], int(columns[2])
                        global_counts[c_id] = global_counts.get(c_id, 0) + reads
                    except ValueError:
                        continue
    return global_counts

def filter_by_global_abundance(current_path, out_path, counts_dict, threshold):
    """Filter current sample based on aggregated global counts."""
    try:
        with open(current_path, 'r', encoding='utf-8') as f_in, \
             open(out_path, 'w', encoding='utf-8', newline='') as f_out:

            reader = csv.reader(f_in, delimiter='\t')
            writer = csv.writer(f_out, delimiter='\t')

            # Transfer Header
            header = next(reader, None)
            if header:
                writer.writerow(header)

            # Apply global filter
            for columns in reader:
                if columns and columns[0] in counts_dict:
                    if counts_dict[columns[0]] >= threshold:
                        writer.writerow(columns)
        return True
    except Exception as e:
        print(f"❌ Error: {e}")
        return False

# --- VARIABLE RETRIEVAL & EXECUTION ---
if __name__ == "__main__":

    path_in = snakemake.input.data
    current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
    path_out = snakemake.output.filtered
    min_total = snakemake.params.reads_threshold

    abundance_dict = get_total_abundance(path_in)

    if filter_by_global_abundance(current_file, path_out, abundance_dict, min_total):
        print(f"✅ Successfully filtered -> {path_out}")
    else:
        raise RuntimeError(f"Filtering failed for {current_file}")
