#!/bin/bash

# 1. Création du dossier pour les logs du cluster
mkdir -p LOGS/snakemake

# 2. Soumission au cluster
# On demande 4 cœurs (-pe thread 4) pour que Snakemake puisse paralléliser
qsub -V -cwd \
     -pe thread 4 \
     -o LOGS/snakemake/ \
     -e LOGS/snakemake/ \
     -q infinit.q \
     -N MicrobExplorer \
     ./run_pipeline.sh "$@"
