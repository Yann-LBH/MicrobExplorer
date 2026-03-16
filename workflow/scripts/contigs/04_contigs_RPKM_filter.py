import csv

# --- SNAKEMAKE CONFIGURATION (VSCode alerts) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(data=["s1_rpkm.tsv", "s2_rpkm.tsv"]), 
        output=SimpleNamespace(rpkm="filtered_rpkm.tsv"), 
        params=SimpleNamespace(rpkm_threshold=0.5), 
        wildcards=SimpleNamespace(sample="s1")
    )

# --- VARIABLE RETRIEVAL ---
path_in_list = snakemake.input.data
# Identify the file corresponding to the current sample wildcard
current_file = [f for f in path_in_list if snakemake.wildcards.sample in f][0]
path_out = snakemake.output.rpkm
threshold = float(snakemake.params.rpkm_threshold)

def get_min_rpkm_across_samples(files):
    """
    Scans all files to find the minimum RPKM value for each contig across the dataset.
    If a contig is missing in one sample, its minimum RPKM is effectively 0.0.
    """
    min_rpkm_map = {}
    is_first_file = True
    
    for path in files:
        current_sample_data = {}
        with open(path, 'r', encoding='utf-8') as f:
            # Using DictReader to access columns by their header names
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                try:
                    contig_id = row['Contig_ID']
                    rpkm_val = float(row['RPKM'])
                    current_sample_data[contig_id] = rpkm_val
                except (ValueError, KeyError):
                    continue

        if is_first_file:
            min_rpkm_map = current_sample_data
            is_first_file = False
        else:
            # Update global map: only keep the intersection of IDs with the minimum value
            for cid in list(min_rpkm_map.keys()):
                if cid in current_sample_data:
                    min_rpkm_map[cid] = min(min_rpkm_map[cid], current_sample_data[cid])
                else:
                    # Contig missing in this sample results in a 0.0 minimum
                    min_rpkm_map[cid] = 0.0
                    
    return min_rpkm_map

def filter_by_min_rpkm(input_path, output_path, min_map, threshold_val):
    """
    Writes the output file keeping only contigs that meet the global minimum RPKM.
    """
    try:
        with open(input_path, 'r', encoding='utf-8') as f_in, \
             open(output_path, 'w', encoding='utf-8', newline='') as f_out:

            reader = csv.DictReader(f_in, delimiter='\t')
            writer = csv.DictWriter(f_out, fieldnames=reader.fieldnames, delimiter='\t')
            
            writer.writeheader()

            for row in reader:
                contig_id = row.get('Contig_ID')
                # Only write row if the global minimum for this ID meets the threshold
                if contig_id and min_map.get(contig_id, 0.0) >= threshold_val:
                    writer.writerow(row)
              
        return True
    except Exception as e:
        print(f"❌ Error during filtering: {e}")
        return False

# --- EXECUTION ---
# 1. Aggregate minimum values across all samples
print(f"📊 Calculating global minimum RPKM across {len(path_in_list)} files...")
global_min_dict = get_min_rpkm_across_samples(path_in_list)

# 2. Filter the current file based on the aggregated data
print(f"📊 Filtering {current_file} (Threshold: {threshold})...")
if filter_by_min_rpkm(current_file, path_out, global_min_dict, threshold):
    print(f"✅ RPKM filtering successful -> {path_out}")
else:
    # Fail the Snakemake job if processing encountered an error
    raise RuntimeError(f"RPKM filtering failed for {current_file}")