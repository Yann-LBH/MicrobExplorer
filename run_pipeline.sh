#!/bin/bash

# Configuration des variables
SNAKEFILE="workflow/snakefile.smk"
CORES=4

echo "========================================="
echo "   Lancement de la pipeline MicrobExplorer"
echo "========================================="

# 1. Activation stricte de Conda pour le script
source $(conda info --base)/etc/profile.d/conda.sh
conda activate snakemake-8.29.0

# 2. Lancement de Snakemake avec les bonnes options
# -p : affiche les commandes shell exécutées
# --use-conda : installe les outils bioinformatiques des règles
snakemake -p -s "$SNAKEFILE" --cores $CORES --use-conda "$@"
