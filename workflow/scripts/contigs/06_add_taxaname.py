import pandas as pd

# --- Compatibility stub (dev only) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(
            data     = "data/s1_rpkm.tsv",
            taxonomy = "data/taxonomy.tsv"
        ),
        output = SimpleNamespace(taxaname="out/s1_annotated.tsv")
    )

def load_taxonomy(tax_path: str) -> pd.Series:
    """
    Charge le fichier taxonomie.
    Format attendu : [0] Class | [1] ID | [2] TaxID | [3] Lineage
    Retourne une Series indexée par Contig_ID.
    """
    df = pd.read_csv(tax_path, sep="\t", header=None,
                     usecols=[1, 3], names=["Contig_ID", "Lineage"])
    df["Lineage"] = df["Lineage"].str.rstrip(";")
    return df.set_index("Contig_ID")["Lineage"]

def run_annotation(input_path: str, output_path: str,
                   taxonomy: pd.Series) -> int:
    """
    Annote les données RPKM avec la taxonomie,
    trie par RPKM décroissant et écrit le résultat.
    Retourne le nombre de contigs annotés.
    """
    df = pd.read_csv(input_path, sep="\t")

    df["Taxonomy"] = df["Contig_ID"].map(taxonomy).fillna("Unclassified")

    df = df.sort_values("RPKM", ascending=False)
    df.to_csv(output_path, sep="\t", index=False)

    return len(df)

# --- Exécution ---
if __name__ == "__main__":
    taxonomy = load_taxonomy(snakemake.input.taxonomy)

    n_kept = run_annotation(
        input_path  = snakemake.input.data,
        output_path = snakemake.output.taxaname,
        taxonomy    = taxonomy
    )
    print(f"✓ {n_kept} contigs annotés -> {snakemake.output.taxaname}")