import pandas as pd

# --- Compatibility stub (dev only) ---
try:
    snakemake
except NameError:
    from types import SimpleNamespace
    snakemake = SimpleNamespace(
        input  = SimpleNamespace(data="data/sample_1.tsv"),
        output = SimpleNamespace(rpkm="results/sample_1_rpkm.tsv"),
    )

def calculate_rpkm(input_path: str, output_path: str) -> int:
    """
    Calcule le RPKM par contig.
    Formule : (Reads * 10^9) / (Length * Total_Mapped_Reads)
    Retourne le nombre de contigs traités.
    """
    df = pd.read_csv(input_path, sep="\t")

    required_cols = {"Reads_Mapped", "Length"}
    missing = required_cols - set(df.columns)
    if missing:
        raise KeyError(f"Colonnes manquantes dans {input_path} : {missing}")

    total_mapped = df["Reads_Mapped"].sum()

    if total_mapped == 0:
        print(f"⚠️  Total mapped reads = 0 pour {input_path}, RPKM mis à NaN.")
        df["RPKM"] = float("nan")
        df.to_csv(output_path, sep="\t", index=False)
        return 0

    df["RPKM"] = ((df["Reads_Mapped"] * 1e9) / (df["Length"] * total_mapped)).round(4)
    df.to_csv(output_path, sep="\t", index=False)
    return len(df)

# --- Exécution ---
if __name__ == "__main__":
    input_file  = snakemake.input.data
    output_file = snakemake.output.rpkm

    n_rows = calculate_rpkm(input_file, output_file)
    print(f"✓ {n_rows} contigs traités -> {output_file}")