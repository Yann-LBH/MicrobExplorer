import pandas as pd

# --- Compatibility stub (dev uniquement) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(all_samples=["data/s1.tsv", "data/s2.tsv"],
                                 current_sample="data/s1.tsv"),
        output = SimpleNamespace(counts="out/s1_filtered.tsv"),
        params = SimpleNamespace(abundance_threshold=10)
    )

# --- Agrégation globale des reads par contig ---
def get_total_abundance(files: list[str]) -> pd.Series:
    """
    Somme les reads par contig_id sur tous les fichiers.
    Retourne une Series indexée par contig_id.
    """
    global_counts = pd.Series(dtype=int)
    for f in files:
        chunk = pd.read_csv(f, sep="\t", usecols=[0, 2])
        global_counts = global_counts.add(
            chunk.groupby(chunk.columns[0])[chunk.columns[2]].sum(),
            fill_value=0
        )
    return global_counts

# --- Filtrage de l'échantillon courant ---
def filter_by_global_abundance(current_path: str, out_path: str,
                                global_counts: pd.Series, threshold: int) -> int:
    """
    Filtre le TSV courant : garde les contigs dont l'abondance globale >= threshold.
    Retourne le nombre de lignes conservées.
    """
    df = pd.read_csv(current_path, sep="\t")
    contig_col = df.columns[0]

    df_filtered = df[df[contig_col].map(global_counts).fillna(0) >= threshold]
    df_filtered.to_csv(out_path, sep="\t", index=False)
    return len(df_filtered)

# --- Exécution ---
if __name__ == "__main__":
    all_files  = snakemake.input.all_samples
    current    = snakemake.input.current_sample
    path_out   = snakemake.output.counts
    threshold  = int(snakemake.params.abundance_threshold)

    global_counts = get_total_abundance(all_files)
    n_kept = filter_by_global_abundance(current, path_out, global_counts, threshold)

    print(f"✓ {n_kept} contigs conservés (seuil abondance={threshold}) -> {path_out}")