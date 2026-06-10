################################################################################
# Project : "MicrobExplorer"
# Script: "Filter contigs by RPKM treshold"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

import pandas as pd

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

def filter_by_min_rpkm(current_path: str, PATH_OUT: str,
                        global_min: pd.Series, RPKM_THRESHOLD: float) -> int:
    """
    Filtre le TSV courant : garde les contigs dont le RPKM min global >= threshold.
    Retourne le nombre de lignes conservées.
    """
    df = pd.read_csv(current_path, sep="\t")
    mask = df["Contig_ID"].map(global_min).fillna(0.0) >= RPKM_THRESHOLD
    df[mask].to_csv(PATH_OUT, sep="\t", index=False)
    return mask.sum()

# --- Exécution ---
if __name__ == "__main__":

    PATH_IN  = snakemake.input.data
    current    = snakemake.input.current_sample
    PATH_OUT   = snakemake.output.rpkm_filtered
    RPKM_THRESHOLD  = float(snakemake.params.rpkm_threshold)

    global_min = get_min_rpkm_across_samples(PATH_IN)
    n_kept     = filter_by_min_rpkm(current, PATH_OUT, global_min, RPKM_THRESHOLD)

    print(f"✓ CONTIGS : Filtered step passed successfully -> {n_kept} contigs kept (threshold={RPKM_THRESHOLD}) -> {PATH_OUT}")