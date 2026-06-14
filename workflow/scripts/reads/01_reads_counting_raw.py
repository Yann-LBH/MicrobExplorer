################################################################################
# Project : "MicrobExplorer"
# Script: "Counting raw reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################


def analyser_kaiju(PATH_IN, PATH_OUT):
    counter = {}

    with open(PATH_IN, "r", encoding="utf-8") as f:
        for line in f:
            # Kaiju sépare généralement les column par des tabulations
            column = line.strip().split("\t")

            # On vérifie si la line est classée ('C') et possède bien un ID taxon
            if len(column) >= 3 and column[0] == "C":
                taxon_id = column[2]  # La 3ème colonne (ex: 35786)

                counter[taxon_id] = counter.get(taxon_id, 0) + 1

    # Écriture des résultats
    if counter:
        with open(PATH_OUT, "w", encoding="utf-8") as f_out:
            # On trie par nombre d'occurrences (du plus grand au plus petit)
            for taxon, total in sorted(
                counter.items(), key=lambda x: x[1], reverse=True
            ):
                f_out.write(f"Taxon {taxon} : {total} reads\n")
        return True
    return False


# ==========================================================================
if __name__ == "__main__":

    PATH_IN = snakemake.input.raw_data
    PATH_OUT = snakemake.output.counted

    success = analyser_kaiju(PATH_IN, PATH_OUT, COUNT_THRESHOLD)
    if success:
        print(f"✓ READS : Counting step passed successfully -> {PATH_OUT}")

    else:
        # Raise an exception so Snakemake knows the script failed
        raise RuntimeError(f"Counting failed for file: {PATH_IN}")
