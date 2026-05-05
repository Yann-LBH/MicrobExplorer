import pandas as pd

# --- Compatibility stub (dev only) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(raw_data="data/sample.tsv"),
        output = SimpleNamespace(counts="out/sample_filtered.tsv"),
        params = SimpleNamespace(length_threshold=1000),
        wildcards = SimpleNamespace(sample="sample")
    )

# --- Filtering ---
def filter_contigs(input_path: str, output_path: str, threshold: int) -> int:
    """
    Reads the input, filters based on `contig_length` >= `threshold`, and writes the output to a TSV file.
    Returns the number of retained contigs.
    """
    df = pd.read_csv(input_path, sep="\t", dtype={"contig_length": int})

    df_filtered = df[df["contig_length"] >= threshold]

    df_filtered.to_csv(output_path, sep="\t", index=False)
    return len(df_filtered)

# --- Execution ---
if __name__ == "__main__":

    path_in   = snakemake.input.raw_data
    path_out  = snakemake.output.counts
    threshold = int(snakemake.params.length_threshold)

    n_kept = filter_contigs(path_in, path_out, threshold)
    print(f"✓ {n_kept} contigs conservés (seuil={threshold}) -> {path_out}")