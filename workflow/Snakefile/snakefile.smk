import pandas as pd
configfile: "../../../config/config.yaml"

metadata = pd.read_csv("config/samples.tsv", sep="\t")

DATA_DIR = config["paths"]["data_raw"]
EXT_READS = config["extensions"]["reads"]
EXT_CONTIGS = config["extensions"]["contigs"]

# Crée une liste vide pour les cibles
final_targets = []

# Ajoute les counts (toujours nécessaires)
final_targets.extend(expand("results/counts/contigs/{sample}_counts.csv", sample=SAMPLES_CONTIGS))

# Taxonomomie
if config.get("run_taxonomy", False):
    final_targets.append("data/taxonomy/formatted_taxo.tsv")

# Contigs processing
if config.get("run_contigs_counts", False):
    final_targets.extend(expand("results/counts/contigs/1.Raw_Counting/{sample}_counts.tsv", sample=SAMPLES_CONTIGS))
if config.get("run_contigs_filtered", False):
    final_targets.extend(expand("results/counts/contigs/2.Filtered/{sample}_filtered.tsv", sample=SAMPLES_CONTIGS))
if config.get("run_contigs_RPKM", False):
    final_targets.extend(expand("results/counts/contigs/3.RPKM/{sample}_RPKM.tsv", sample=SAMPLES_CONTIGS))
if config.get("run_contigs_RPKM_filter", False):
    final_targets.extend(expand("results/counts/contigs/4.RPKM_Filtered/{sample}_RPKM_filtered.tsv", sample=SAMPLES_CONTIGS))
if config.get("run_contigs_intersec", False):
    final_targets.extend(expand("results/counts/contigs/5.Intersection/{sample}_counts.tsv", sample=SAMPLES_CONTIGS))

# Plots
if config.get("run_barplot", False):
    final_targets.append("results/barplot.png")
if config.get("run_heatmap", False):
    final_targets.append("results/heatmap.png")
if config.get("run_pca", False):
    final_targets.append("results/pca_plot.png")

if not SAMPLES_READS and not SAMPLES_CONTIGS:
    raise ValueError(f"Aucun fichier trouvé dans {DATA_DIR} avec les extensions {EXT_READS} ou {EXT_CONTIGS}")

# ---------------------------------------------------------------------------- Final targets
rule all:
    input: final_targets

# ---------------------------------------------------------------------------- Taxanomie table processing        
rule download_ncbi:
    output:
        temp("data/taxonomy/ncbi_dump.tar.gz") # Sera supprimé après la règle suivante
    shell:
        "wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O {output}"

rule convert_taxo:
    input:
        "data/taxonomy/ncbi_dump.tar.gz"
    output:
        "data/taxonomy/formatted_taxo.tsv" # Ce fichier-là est conservé
    script:
        "scripts/convert_ncbi.py"

# ---------------------------------------------------------------------------- Processing
# ----------------- Step 1
rule contigs_counts:
    input:
        raw_data = "data/raw/{sample}.tsv"
    output:
        counts = "results/counts/contigs/1.Raw_Counting/{sample}_counts.tsv"
        step_checking = "report/check_steps/{sample}_abundance_check.tsv"
    log:
        "results/logs/contigs/{sample}_filter.log"
    params : 
        length_threshold = config["contigs"]["length_threshold"]
    conda:
        "envs/python_env.yaml"
    script:
        "../scripts/contigs/01_contigs_counting_raw.py"
# ----------------- Step 2
rule contigs_filtered:
    input:
        data = expand("results/counts/contigs/1.Raw_Counting/{sample}_counts.tsv", sample=SAMPLES_CONTIGS)
    output:
        filtered ="results/counts/contigs/2.Filtered/{sample}_filtered.tsv"
        step_checking = "report/check_steps/{sample}_abundance_check.tsv"
    log:
        "results/logs/contigs/{sample}_filtered.log"
    params : 
        reads_threshold = config["contigs"]["reads_threshold"]
    conda:
        "envs/python_env.yaml"
    script:
        "../scripts/contigs/02_contigs_filter.py"
# ----------------- Step 3
rule contigs_RPKM:
    input:
        data = "results/counts/contigs/2.Filtered/{sample}_filtered.tsv"
    output:
        rpkm ="results/counts/contigs/3.RPKM/{sample}_RPKM.tsv"
        step_checking = "report/check_steps/{sample}_abundance_check.tsv"
    log:
        "results/logs/contigs/{sample}_RPKM.log"
    conda:
        "envs/python_env.yaml"
    script:
        "../scripts/contigs/03_contigs_RPKM.py"
# ----------------- Step 4
rule contigs_RPKM_filter:
    input:
        data = expand("results/counts/contigs/3.RPKM/{sample}_RPKM.tsv", sample=SAMPLES_CONTIGS)
    output:
        rpkm_filtered ="results/counts/contigs/4.RPKM_Filtered/{sample}_RPKM_filtered.tsv"
        step_checking = "report/check_steps/{sample}_abundance_check.tsv"
    log:
        "results/logs/contigs/{sample}_RPKM_trimming.log"
    params : 
        rpkm_threshold = config["contigs"]["rpkm_threshold"]
    conda:
        "envs/python_env.yaml"
    script:
        "../scripts/contigs/04_contigs_RPKM_filter.py"
# ----------------- Step 5
rule contigs_intersec:
    input:
        "results/counts/contigs/4.RPKM_Filtered/{sample}_RPKM_filtered.tsv"
    output:
        "results/counts/contigs/5.Intersection/{sample}_counts.tsv"
    log:
        "results/logs/contigs/{sample}_intersection.log"
    conda:
        "envs/python_env.yaml"
    script:
        "../scripts/contigs/05_contigs_intersec_filtered.py"
# ---------------------------------------------------------------------------- Plotting
rule barplot:
    input:
        counts = "results/counts.csv"
    output:
        barplot = "results/barplot.pdf"
    params:
        barplot_color = config["plots"]["barplot"]["color_by"]
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/08_barplot_DNA.R"

rule plot_heatmap:
    input:
        data = "results/counts.csv"
    output:
        heatmap = "results/heatmap.png"
    params:
        color = config["plots"]["heatmap"]["color_by"]
        size = config["plots"]["heatmap"]["size"]
        method = config["plots"]["heatmap"]["method"]
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/09_heatmap.R"

rule plot_pca:
    input:
        counts = "results/counts.csv"
    output:
        pca = "results/pca_plot.png"
    params:
        pca_color = config["plots"]["pca"]["color_by"]
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/10_pca.R"    
        
        # output: directory("results/taxonomie_complete") # On déclare le dossier entier comme sortie
        
rule process_reads:
    input: "data/raw/reads.kaiju"
    output: "results/intermediaire/reads_taxo_clean.csv"
    script: "scripts/01_process_reads_taxo.py"

rule process_contigs:
    input: "data/raw/contigs.kaiju"
    output: "results/intermediaire/contigs_taxo_clean.csv"
    script: "scripts/01_process_contigs_taxo.py"

rule merge_data:
    input:
        r = "results/intermediaire/reads_taxo_clean.csv",
        c = "results/intermediaire/contigs_taxo_clean.csv"
    output: "results/final_counts.csv"
    script: "scripts/02_merge_all_sources.py"

# ---------------------------------------------------------------------------- Quality Check
rule qc:
    input: ######
    ouput: "report/qc_contigs.tsv"
    script: "scripts/utils/B_contigsqc.py"
