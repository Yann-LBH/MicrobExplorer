# ==========================================================================
# Pipeline Metagenomics - Configuration
# ==========================================================================

import pandas as pd
configfile: "../../../config/config.yaml"

# Lecture des échantillons
metadata = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = metadata['sample'].tolist() 

# --- Configuration des cibles (Targets) ---
final_targets = []

# Taxonomie
if config.get("run_taxonomy"):
    final_targets.append("data/taxonomy/formatted_taxo.tsv")

# Contigs processing
# Dictionnaire des étapes de processing
CONTIGS_STEPS = {
    "1.Raw_Counting":   ["01_contigs_counting_raw.py",  "data/raw/{sample}.tsv"],
    "2.Filtered":       ["02_contigs_filter.py",        "results/treatment/contigs/1.Raw_Counting/{sample}_counts.tsv"],
    "3.rpkm":           ["03_contigs_rpkm.py",          "results/treatment/contigs/2.Filtered/{sample}_filtered.tsv"],
    "4.rpkm_Filtered":  ["04_contigs_rpkm_filter.py",   "results/treatment/contigs/3.rpkm/{sample}_rpkm.tsv"],
    "5.Intersection":   ["05_contigs_intersec.py",      "results/treatment/contigs/4.rpkm_Filtered/{sample}_rpkm_filtered.tsv"],
    "6.Final_Results":  ["06_contigs_add_taxaname.py",  "results/treatment/contigs/5.Intersection/{sample}_intersec.tsv"]
}

READS_STEPS = {
    "1.Raw_Counting":   ["01_reads_counting_raw.py",  "data/raw/{sample}.tsv"],
    "2.Filtered":       ["02_reads_filter.py",        "results/treatment/reads/1.Raw_Counting/{sample}_counts.tsv"],
    "3.Final_Results":  ["03_reads_add_taxaname.py",  "results/treatment/reads/2.Filtered/{sample}_filtered.tsv"]
}

# Plots
PLOT_CONFIG = {
    "barplot": ["07_barplot.R", "data/processed/barplots_data.rds"],
    "stacked_barplot": ["07_barplot.R", "data/processed/stacked_barplot_data.rds"],
    "heatmap": ["08_heatmap.R", "data/processed/all_samples_consolidated.rds"],
    "pca":     ["09_ACP.R",     "data/processed/all_samples_consolidated.rds"],
    "volcano": ["10_volcano.R", "data/processed/deseq2_results.rds"],
}

PLOT_CONFIG = {
    "permanova": ["07_barplot.R", "data/processed/barplots_data.rds"],
    "deseq2": ["07_barplot.R", "data/processed/barplots_data.rds"]
}

PLOT_TYPES = ["barplot", "heatmap", "pca", "volcano", "permanova"]
for plot_type in PLOT_TYPES:
    if config.get(f"run_{plot_type}"):
        final_targets.append(f"results/plots/{plot_type}.png")

if not SAMPLES_READS and not SAMPLES_CONTIGS:
    raise ValueError(f"Aucun fichier trouvé dans {DATA_DIR} avec les extensions {EXT_READS} ou {EXT_CONTIGS}")

for opt, path in CONTIGS_STEPS.items():
    if config.get(opt): final_targets.extend(expand(f"results/counts/contigs/{path}", sample=SAMPLES))

# ---------------------------------------------------------------------------- Final targets
rule all:
    input: final_targets

# ---------------------------------------------------------------------------- Taxanomie table processing        
rule taxonomy_full:
    output: "data/taxonomy/formatted_taxo.tsv"
    shell: 
        "wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O - | tar -xz && "
        "python scripts/convert_ncbi.py"

# ---------------------------------------------------------------------------- Processing
rule contigs_processing:
    input: lambda w: CONTIGS_STEPS[w.step][1].format(sample=w.sample)
    output: 
        result = "results/counts/contigs/{step}/{sample}_{name}.tsv",
    params: threshold = config["contigs"]["length_threshold"]
    conda: "envs/python_env.yaml"
    script: lambda w: f"../scripts/contigs/{CONTIGS_STEPS[w.step][0]}"

rule reads_processing:
    input: lambda w: RREADS_STEPS[w.step][1].format(sample=w.sample)
    output: 
        counts = "results/counts/contigs/{step}/{sample}_{name}.tsv",
    params: threshold = config["contigs"]["length_threshold"]
    conda: "envs/python_env.yaml"
    script: lambda w: f"../scripts/contigs/{RREADS_STEPS[w.step][0]}"
# ---------------------------------------------------------------------------- Plotting
# --- 1. Préparation pour les Barplots (données Step 2) ---
rule prepare_barplot_data:
    input:  expand("results/counts/contigs/2.Filtered/{sample}_filtered.tsv", sample=SAMPLES)
    output: "data/processed/barplots_data.rds"
    script: "../scripts/utils/prep_barplots.R"

# --- 2. Préparation pour Analyse Différentielle (Volcano) ---
rule run_deseq2:
    input:  expand("results/counts/contigs/6.Final_Results/{sample}_intersec.tsv", sample=SAMPLES)
    output: "data/processed/deseq2_results.rds"
    script: "../scripts/utils/run_deseq2.R"

# Règle générique pour les plots
rule generate_plots:
    input:
        rds = lambda w: PLOT_CONFIG[w.plot_name][1]
    output:
        fig = "results/plots/{plot_name}.png"
    params: 
        color   = lambda w: config["plots"].get(w.plot_name, {}).get("color_by"),
        size    = lambda w: config["plots"].get(w.plot_name, {}).get("size"),
        method  = lambda w: config["plots"].get(w.plot_name, {}).get("method")
    conda:
        "envs/r_env.yaml"
    script:
        lambda w: f"../scripts/plotting/{PLOT_CONFIG[w.plot_name][0]}"
        "Rscript -e 'renv::restore(); source(\"scripts/plotting/07_barplot.R\")'"


rule all:
    input:
        "results/plots_parameters/"

rule generate_physico_plots:
    input:
        "data/parameters.xlsx"
    output:
        directory("results/plots_parameters/")
    script:
        "scripts/plot_parameters.R"








rule generate_plots:
    input:  lambda w: PLOTS[w.plot_name][1]
    output: fig = "results/plots/{plot_name}.png"
    params: color   = lambda w: config["plots"].get(w.plot_name, {}).get("color_by"),
            size    = lambda w: config["plots"].get(w.plot_name, {}).get("size"),
            method  = lambda w: config["plots"].get(w.plot_name, {}).get("method")
    conda:  "envs/r_env.yaml"
    script: "scripts/plotting/{w.plot_name}.R"
        



        
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


# Dans ton Snakefile fin traitement avant script R
rule aggregate_filtered_contigs:
    input:
        expand("out/{sample}_filtered.tsv", sample=SAMPLES)
    output:
        "Data/Parquet/contigs_filtered.parquet"
    run:
        import pandas as pd

        pd.concat(
            [pd.read_csv(f, sep="\t").assign(sample=f.split("/")[-1].replace("_filtered.tsv",""))
             for f in input],
            ignore_index=True
        ).to_parquet(output[0], index=False)

# Script 2 contigs
rule filter_abundance:
    input:
        all_samples   = expand("out/{sample}_step1.tsv", sample=SAMPLES),
        current_sample = "out/{sample}_step1.tsv"
    output:
        counts = "out/{sample}_step2.tsv"
    params:
        abundance_threshold = 10
    script:
        "scripts/filter_abundance.py"

rule pca:
    input:
        parquet  = expand("Data/Parquet/{sample}_annotated.parquet", sample=SAMPLES),
        metadata = "metadata/metadata.xlsx",
        physico  = "Physico-chimique/Thermophiles.xlsx"
    output:
        pdf     = "Graphique/ACP/ACP_results.pdf",
        parquet = "Data/Parquet/pca_results.parquet",
        csv     = "Graphique/ACP/contributions.csv"
    params:
        dim_x        = 1,
        dim_y        = 2,
        physico_cols = ["pH", "AGV mg/l"],
        n_top        = 10
    script:
        "scripts/pca.R"