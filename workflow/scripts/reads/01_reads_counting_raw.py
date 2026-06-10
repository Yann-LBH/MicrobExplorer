################################################################################
# Project : "MicrobExplorer"
# Script: "Counting raw reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

def analyser_kaiju(path_in, path_out):
    counter = {}
    
    with open(path_in, 'r', encoding='utf-8') as f:
        for line in f:
            # Kaiju sépare généralement les column par des tabulations
            column = line.strip().split('\t')
            
            # On vérifie si la line est classée ('C') et possède bien un ID taxon
            if len(column) >= 3 and column[0] == 'C':
                taxon_id = column[2] # La 3ème colonne (ex: 35786)
                
                counter[taxon_id] = counter.get(taxon_id, 0) + 1

    # Écriture des résultats
    if counter:
        with open(path_out, 'w', encoding='utf-8') as f_out:
            # On trie par nombre d'occurrences (du plus grand au plus petit)
            for taxon, total in sorted(counter.items(), key=lambda x: x[1], reverse=True):
                f_out.write(f"Taxon {taxon} : {total} reads\n")
        return True
    return False

# ==========================================================================
if __name__ == "__main__":

    path_in     = snakemake.input.raw_data
    path_out    = snakemake.output.counted
    current     = snakemake.input.current_sample

    print(f"✓ Counting step passed -> {path_out}")