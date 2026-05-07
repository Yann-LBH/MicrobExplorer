################################################################################
# Project : "MicrobExplorer"
# Script: "Filter for contigs"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

# --- Agrégation globale des reads par contig ---
def get_total_abundance(path_in: list[str]) -> pd.Series:
    """
    Somme les reads par contig_id sur tous les fichiers.
    Retourne une Series indexée par contig_id.
    """
    global_counts = pd.Series(dtype=int)
    for f in path_in:
        chunk = pd.read_csv(f, sep="\t", usecols=[0, 2])
        global_counts = global_counts.add(
            chunk.groupby(chunk.columns[0])[chunk.columns[2]].sum(), # peut etre 0 et 1 ici 
            fill_value=0
        )
    return global_counts

# --- Filtrage de l'échantillon courant ---
def filter_by_global_abundance(current: str, path_out: str,
                                global_counts: pd.Series, threshold: int) -> int:
    """
    Filtre le TSV courant : garde les contigs dont l'abondance globale >= threshold.
    Retourne le nombre de lignes conservées.
    """
    df = pd.read_csv(current, sep="\t")
    contig_col = df.columns[0]

    df_filtered = df[df[contig_col].map(global_counts).fillna(0) >= threshold]
    df_filtered.to_csv(path_out, sep="\t", index=False)
    return len(df_filtered)

# --- Exécution ---
if __name__ == "__main__":
    path_in    = snakemake.input.data
    current    = snakemake.input.current_sample
    path_out   = snakemake.output.filtered
    threshold  = int(snakemake.params.abundance_threshold)

    global_counts = get_total_abundance(path_in)
    n_kept = filter_by_global_abundance(current, path_out, global_counts, threshold)

    print(f"✓ {n_kept} contigs kept (threshold abundance={threshold}) -> {path_out}")