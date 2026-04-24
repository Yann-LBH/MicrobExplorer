import pandas as pd

# --- Compatibility stub (dev only) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(data=["data/s1.tsv", "data/s2.tsv"],
                                 current_sample="data/s1.tsv"),
        output = SimpleNamespace(rpkm="out/s1_filtered.tsv"),
        params = SimpleNamespace(rpkm_threshold=0.1)
    )

def get_min_rpkm_across_samples(files: list[str]) -> pd.Series:
    """
    Retourne le RPKM minimum par contig sur tous les échantillons.
    Un contig absent d'un échantillon vaut 0.0 pour cet échantillon.
    """
    frames = [
        pd.read_csv(f, sep="\t", usecols=["Contig_ID", "RPKM"])
          .set_index("Contig_ID")["RPKM"]
        for f in files
    ]
    # concat en colonnes puis min par ligne — gère les NaN (contig absent) comme 0.0
    return pd.concat(frames, axis=1).fillna(0.0).min(axis=1)

def filter_by_min_rpkm(current_path: str, output_path: str,
                        global_min: pd.Series, threshold: float) -> int:
    """
    Filtre le TSV courant : garde les contigs dont le RPKM min global >= threshold.
    Retourne le nombre de lignes conservées.
    """
    df = pd.read_csv(current_path, sep="\t")
    mask = df["Contig_ID"].map(global_min).fillna(0.0) >= threshold
    df[mask].to_csv(output_path, sep="\t", index=False)
    return mask.sum()

# --- Exécution ---
if __name__ == "__main__":
    all_files  = snakemake.input.data
    current    = snakemake.input.current_sample
    path_out   = snakemake.output.rpkm
    threshold  = float(snakemake.params.rpkm_threshold)

    global_min = get_min_rpkm_across_samples(all_files)
    n_kept     = filter_by_min_rpkm(current, path_out, global_min, threshold)

    print(f"✓ {n_kept} contigs conservés (seuil RPKM={threshold}) -> {path_out}")