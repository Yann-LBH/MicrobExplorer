################################################################################
# Project : "MicrobExplorer"
# Script  : "Rank abundance from Kaiju reads"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
################################################################################

import os
import sys
import shutil
import zipfile
import argparse
from collections import defaultdict

# ---------------------------------------------------------------------------
# Rank column index in rankedlineage.dmp (pipe-separated)
# field 0  = tax_id
# field 1  = scientific_name
# field 2  = species        field 3 = genus
# field 4  = family         field 5 = order
# field 6  = class          field 7 = phylum
# field 8  = kingdom        field 9  = domain
# ---------------------------------------------------------------------------
RANK_COLUMN = {
    "domain":  9,
    "kingdom": 8,
    "phylum":  7,
    "class":   6,
    "order":   5,
    "family":  4,
    "genus":   3,
    "species": 2,
}

TOP_N = 19  # Number of top taxa to keep before collapsing the rest into "Others"


def build_tax_map(zip_path: str, dmp_name: str, rank_col: int) -> dict[str, str]:
    """
    Extract tax_id → rank_name mapping from rankedlineage.dmp inside the NCBI zip.

    Parameters
    ----------
    zip_path  : path to new_taxdump.zip
    dmp_name  : name of the dmp file inside the zip
    rank_col  : column index of the desired rank (see RANK_COLUMN)

    Returns
    -------
    dict mapping tax_id (str) → rank label (str)
    """
    rank_label = next(k for k, v in RANK_COLUMN.items() if v == rank_col)
    unknown_label = f"Unknown {rank_label.capitalize()}"

    tax_map: dict[str, str] = {}

    with zipfile.ZipFile(zip_path) as z, z.open(dmp_name) as fh:
        for raw in fh:
            parts = raw.decode("utf-8").split("|")
            if len(parts) > rank_col:
                tax_id = parts[0].strip()
                name   = parts[rank_col].strip()
                tax_map[tax_id] = name if name else unknown_label

    return tax_map


def count_reads(kaiju_path: str, tax_map: dict[str, str]) -> tuple[dict[str, int], int]:
    """
    Count reads per rank label from a single Kaiju output file.

    Parameters
    ----------
    kaiju_path : path to a .kaijuNR file
    tax_map    : tax_id → rank_label mapping

    Returns
    -------
    (counter, total_reads)
    """
    counter: dict[str, int] = defaultdict(int)
    total = 0

    with open(kaiju_path, encoding="utf-8") as fh:
        for line in fh:
            total += 1
            cols = line.rstrip("\n").split("\t")
            if cols[0] == "C":
                label = tax_map.get(cols[2], "Unknown/Other")
            else:
                label = "Unclassified"
            counter[label] += 1

    return dict(counter), total


def write_output(out_path: str, counter: dict[str, int], total: int,
                 rank_label: str) -> None:
    """
    Write abundance TSV: top TOP_N taxa + one aggregated "Others" row.

    Parameters
    ----------
    out_path   : destination file path
    counter    : {rank_label: read_count}
    total      : total read count (classified + unclassified)
    rank_label : column header name for the rank (e.g. "Phylum")
    """
    sorted_items = sorted(counter.items(), key=lambda x: x[1], reverse=True)
    top          = sorted_items[:TOP_N]
    others_sum   = sum(c for _, c in sorted_items[TOP_N:])

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(f"{rank_label};Reads;Abundance_Percent\n")
        for name, count in top:
            fh.write(f"{name};{count};{count / total * 100:.2f}%\n")
        if others_sum:
            fh.write(f"Others;{others_sum};{others_sum / total * 100:.2f}%\n")


def process_folder(zip_path: str, dmp_name: str, data_folder: str,
                   out_folder: str, rank: str) -> None:
    """
    Run the full abundance pipeline for every .kaijuNR file in data_folder.

    Parameters
    ----------
    zip_path    : path to new_taxdump.zip
    dmp_name    : dmp filename inside the zip
    data_folder : folder containing .kaijuNR files
    out_folder  : output folder for abundance TSVs
    rank        : rank name (key of RANK_COLUMN)
    """
    rank = rank.lower()
    if rank not in RANK_COLUMN:
        sys.exit(f"Unknown rank '{rank}'. Valid options: {list(RANK_COLUMN)}")

    rank_col = RANK_COLUMN[rank]

    # Clean and recreate output folder
    if os.path.exists(out_folder):
        shutil.rmtree(out_folder)
    os.makedirs(out_folder)

    print(f"Building taxonomy map for rank '{rank}' (column {rank_col})…")
    tax_map = build_tax_map(zip_path, dmp_name, rank_col)
    print(f"  {len(tax_map):,} tax IDs loaded.")

    kaiju_files = [f for f in os.listdir(data_folder) if f.endswith(".kaijuNR")]
    if not kaiju_files:
        print(f"No .kaijuNR files found in '{data_folder}'.", file=sys.stderr)
        return

    for filename in kaiju_files:
        in_path  = os.path.join(data_folder, filename)
        out_name = f"abundance_{rank}_{filename.replace('.kaijuNR', '.tsv')}"
        out_path = os.path.join(out_folder, out_name)

        print(f"Processing {filename}…")
        counter, total = count_reads(in_path, tax_map)

        if total == 0:
            print(f"  ⚠ No reads found in {filename}.")
            continue

        write_output(out_path, counter, total, rank.capitalize())
        print(f"  ✓ {total:,} reads → {out_path}")

    print("Done.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compute read abundance at any NCBI taxonomic rank."
    )
    p.add_argument("--zip",    default="A.Taxa_Table_NCBI/new_taxdump.zip",
                   help="Path to new_taxdump.zip")
    p.add_argument("--dmp",    default="rankedlineage.dmp",
                   help="Name of the dmp file inside the zip")
    p.add_argument("--data",   default="Data",
                   help="Folder containing .kaijuNR files")
    p.add_argument("--out",    default=None,
                   help="Output folder (default: Abundance/<rank>)")
    p.add_argument("--rank",   default="phylum",
                   choices=list(RANK_COLUMN),
                   help="Taxonomic rank to aggregate on")
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    out_folder = args.out or os.path.join("Abundance", args.rank.capitalize())
    process_folder(args.zip, args.dmp, args.data, out_folder, args.rank)