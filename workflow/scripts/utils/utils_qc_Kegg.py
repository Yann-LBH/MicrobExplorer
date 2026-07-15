################################################################################
# Project : "MicrobExplorer"
# Script: "Utils for Kegg QC"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import os

def get_kegg_ids_from_gff3(path):
    """Extract unique KEGG (KO) IDs from a standard GFF3 file attributes column."""
    kegg_ids = set()
    try:
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                # Skip comments and empty lines
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) < 9:
                    continue
                
                attributes = parts[8]
                # Look for KEGG or KO annotations in the attributes column (9th column)
                if "KEGG:" in attributes or "Dbxref=KEGG:" in attributes:
                    # Common GFF3 parsing logic for cross-references
                    for attr in attributes.split(';'):
                        if "KEGG:" in attr:
                            ko_id = attr.split(':', 1)[1].strip()
                            # Ensure it looks like a KO number (e.g., K00844)
                            if ko_id.startswith('K'):
                                kegg_ids.add(ko_id.split(',')[0]) # Handle potential list
    except Exception as e:
        print(f"Error reading GFF3 file {path}: {e}")
    return list(kegg_ids)

def get_kegg_ids_from_tsv(path):
    """Read a standard sample-specific TSV file and return a list of unique KEGG IDs."""
    kegg_ids = set()
    try:
        with open(path, 'r', encoding='utf-8') as f:
            next(f, None) # Skip header
            for line in f:
                if line.strip():
                    parts = line.split('\t')
                    if parts[0].strip():
                        kegg_ids.add(parts[0].strip())
    except Exception as e:
        print(f"Error reading TSV file {path}: {e}")
    return list(kegg_ids)

def run_full_qc(snakemake, steps_config):
    """Main execution function for KEGG QC handling a single raw GFF3 and sample TSVs."""
    DATA = snakemake.input.data
    qc_results = {}
    
    # Get the list of all samples from the snakemake object params or globals
    # Assuming snakemake.params.samples or we can infer them during cleanup
    all_samples = set()

    # First pass: identify all sample names from the TSV files present in the input
    for path in DATA:
        filename = os.path.basename(path)
        if not filename.endswith(".gff3"):
            sample_name = filename
            for item in [
                "extracted_",
                "intersected_",
                "deseq2_",
                "standardized_",
                "stand_aggreg_",
                "annotated_",
                "_kegg"
                ".tsv", 
            ]:
                sample_name = sample_name.replace(item, "")
            all_samples.add(sample_name)

    # Second pass: Process the files
    for path in DATA:
        filename = os.path.basename(path)
        normalized_path = os.path.normpath(path)

        # 1. Step Detection
        step_key = "Unknown"
        for key, (folder, pattern) in steps_config.items():
            normalized_folder = os.path.normpath(folder)
            if normalized_folder in normalized_path:
                step_key = key
                break

        # 2. Extract Counts based on file type
        if filename.endswith(".gff3"):
            # It's the global raw co-assembly file
            total_raw_kegg = len(get_kegg_ids_from_gff3(path))
            
            # Assign this identical baseline count to ALL known samples
            for sample_name in all_samples:
                if sample_name not in qc_results:
                    qc_results[sample_name] = {}
                qc_results[sample_name][step_key] = total_raw_kegg
        else:
            # It's a standard sample-specific TSV file
            sample_name = filename
            for item in [
                "extracted_",
                "intersected_",
                "deseq2_",
                "standardized_",
                "stand_aggreg_",
                "annotated_",
                "_kegg",
                ".tsv",
            ]:
                sample_name = sample_name.replace(item, "")
                
            ids = get_kegg_ids_from_tsv(path)
            
            if sample_name not in qc_results:
                qc_results[sample_name] = {}
            qc_results[sample_name][step_key] = len(ids)

    return qc_results