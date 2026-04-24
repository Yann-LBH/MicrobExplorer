import pandas as pd

# --- Compatibility stub (dev only) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(
            abund    = "data/s1_abund_filtered.tsv",
            rpkm     = "data/s1_rpkm_filtered.tsv",
            source   = "data/s1_rpkm.tsv"
        ),
        output = SimpleNamespace(out="out/s1_final.tsv")
    )

def get_union_and_extract(path_abund: str, path_rpkm_filt: str,
                           path_source: str, path_out: str) -> int:
    """
    Union des Contig_ID des deux fichiers filtrés,
    puis extraction des lignes correspondantes depuis le fichier source.
    Retourne le nombre de contigs extraits.
    """
    ids_abund = pd.read_csv(path_abund,    sep="\t", usecols=["Contig_ID"])["Contig_ID"]
    ids_rpkm  = pd.read_csv(path_rpkm_filt, sep="\t", usecols=["Contig_ID"])["Contig_ID"]

    target_ids = set(ids_abund).union(ids_rpkm)

    df_source  = pd.read_csv(path_source, sep="\t")
    df_out     = df_source[df_source["Contig_ID"].isin(target_ids)]
    df_out.to_csv(path_out, sep="\t", index=False)

    return len(df_out)

# --- Exécution ---
if __name__ == "__main__":
    n_kept = get_union_and_extract(
        path_abund    = snakemake.input.abund,
        path_rpkm_filt = snakemake.input.rpkm,
        path_source   = snakemake.input.source,
        path_out      = snakemake.output.out
    )
    print(f"✓ {n_kept} contigs extraits -> {snakemake.output.out}")