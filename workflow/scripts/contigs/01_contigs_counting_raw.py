import csv

# --- SNAKEMAKE CONFIGURATION (VSCode alerts) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(), 
        output=SimpleNamespace(), 
        params=SimpleNamespace(), 
        wildcards=SimpleNamespace()
    )

# --- VARIABLE RETRIEVAL ---
path_in = snakemake.input.raw_data
current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
path_out = snakemake.output.counts
min_size = int(snakemake.params.length_threshold)
sample_id = snakemake.wildcards.sample

def filter_tsv_contigs(path_in, path_out, threshold):
    """
    Filters a TSV file of contigs based on a minimum length threshold.
    """
    kept_rows = []
    
    try:
        with open(path_in, 'r', encoding='utf-8') as f_in:
            # TSV reading in dictionary mode
            reader = csv.DictReader(f_in, delimiter='\t')
            fieldnames = reader.fieldnames
            
            for row in reader:
                try:
                    # Checking 'Length' column
                    if int(row['Length']) >= threshold:
                        kept_rows.append(row)
                except (ValueError, KeyError):
                    # Skip the row if 'Length' is missing or not an integer
                    continue

        # --- WRITING THE FILTERED FILE ---
        with open(path_out, 'w', encoding='utf-8', newline='') as f_out:
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            writer.writerows(kept_rows)

        return True

    except Exception as e:
        print(f"❌ Error during the treatment of {path_in} : {e}")
        return False

# --- EXECUTION ---
print(f"📊 Analyze of {path_in} (Threshold: {min_size})...")

if filter_tsv_contigs(current_file, path_in, path_out, min_size):
    print(f"✅ Filtered successfully -> {path_out}")
else:
    # ERROR HANDLING 
    raise RuntimeError(f"❌ Failed to count contigs for {current_file}")
