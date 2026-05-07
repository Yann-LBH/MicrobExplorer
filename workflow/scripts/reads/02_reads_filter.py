################################################################################
# Project : "MicrobExplorer"
# Script: "Filter for reads"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

def trimming_kaiju(path_in, path_out):
    
    try:
        with open(path_in, 'r', encoding='utf-8') as f_in, \
             open(path_out, 'w', encoding='utf-8') as f_out:

            #garde les reads >= 100
            for line in f_in:
            # On découpe la ligne par espaces
                columns = line.strip().split(' ')
            
                if len(columns) >= 4:
                    try:
                        # Conversion du nombre de reads en entier
                        number_reads = int(columns[3]) 
                    
                        # On ne garde que si c'est >= 100
                        if number_reads >= 100:
                            f_out.write(line)
                    except ValueError:
                        # En cas d'en-tête (texte au lieu de nombre), on ignore l'erreur
                        continue
            #garde les 100 premières lignes seulement ( les plus abondants)
            #for i, lines in enumerate(f_in):
                #if i >= 100:
                    #break
                #f_out.write(lines)

        return True
    except Exception as e:
            print(f"❌ Erreur sur le fichier {path_in}: {e}")
            return False
    
# ==========================================================================
if __name__ == "__main__":

    path_in        = snakemake.input.data
    path_out     = snakemake.output.filtered

    print(f"✓ Filtering step passed -> {path_out}")
    