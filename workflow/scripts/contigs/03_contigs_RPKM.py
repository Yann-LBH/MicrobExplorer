import csv

# --- SNAKEMAKE CONFIGURATION (VSCode alerts) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input=SimpleNamespace(), 
        output=SimpleNamespace(),
        wildcards=SimpleNamespace()
    )

# --- VARIABLE RETRIEVAL ---
path_in = snakemake.input.data
current_file = [f for f in path_in if snakemake.wildcards.sample in f][0]
path_out = snakemake.output.rpkm

def calculate_rpkm_for_file(p_in, p_out):
    """
    Calculates RPKM for each contig in the file.
    Formula: (Mapped_Reads * 10^9) / (Contig_Length * Total_Mapped_Reads_in_Sample)
    """
    try:
        # Step 1: Calculate Total Mapped Reads (N) for the normalization factor
        total_mapped_sample = 0
        with open(p_in, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                try:
                    total_mapped_sample += int(row['Reads_Mapped'])
                except (ValueError, KeyError):
                    continue
        
        # Avoid division by zero if the sample has no reads
        if total_mapped_sample == 0:
            print("⚠️ Warning: Total mapped reads is zero. Cannot calculate RPKM.")
            return False

        # Step 2: Calculate RPKM and write output
        with open(p_in, 'r', encoding='utf-8') as f_in, \
             open(p_out, 'w', encoding='utf-8', newline='') as f_out:
            
            reader = csv.DictReader(f_in, delimiter='\t')
            # Prepare output header: original columns + RPKM
            fieldnames = reader.fieldnames + ['RPKM']
            writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()

            for row in reader:
                try:
                    reads = int(row['Reads_Mapped'])
                    length = int(row['Length'])
                    
                    # Apply RPKM Formula
                    rpkm_val = (reads * 10**9) / (length * total_mapped_sample)
                    row['RPKM'] = f"{rpkm_val:.4f}"
                    writer.writerow(row)
                except (ValueError, KeyError, ZeroDivisionError):
                    # Skip lines with invalid numeric data or zero length
                    continue

        return True

    except Exception as e:
        print(f"❌ Error: {e}")
        return False

# --- EXECUTION ---
print(f"📊 Calculating RPKM for: {current_file}")
if calculate_rpkm_for_file(current_file, path_out):
    print(f"✅ RPKM calculation successful -> {path_out}")
else:
    # Raise error so Snakemake stops the pipeline
    raise RuntimeError(f"RPKM calculation failed for {current_file}")