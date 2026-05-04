import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# ==========================================================================
# Échantillons
# ==========================================================================

data           = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES        = data["sample"].tolist()
CONTIGS_FILES  = dict(zip(data["sample"], data["contigs_file"]))
READS_FILES    = dict(zip(data["sample"], data["reads_file"]))
KEGG_FILES     = dict(zip(data["sample"], data["kegg_file"]))

# ==========================================================================
# Cibles finales
# ==========================================================================
def get_targets():
    targets = []

    # --- Taxonomy ---
    if config["run_taxonomy"]:
        targets.append("data/taxonomy/mapping_taxons.csv") #####

    # --- Contigs ---
    if config["run_contigs"]:
        targets.extend(expand(
            "results/treatment/contigs/6.Final_Results/{sample}_annotated.tsv",
            sample=SAMPLES
        ))

    # --- Reads ---
    if config["run_reads"]:
        targets.extend(expand(
            "results/treatment/reads/2.Final_Results/{sample}_annotated.tsv",
            sample=SAMPLES
        ))

    # --- QC ---
    if config["run_qc"]:
        targets.append("results/qc/Rapport_QC_Final.pdf")

    # --- Plots ---
    if config["run_plots"]:
        targets += [
            "results/plots/heatmap/Heatmaps.pdf",
            "results/plots/barplot/Barplots.pdf",
            source=["contigs", "reads"]
            "results/plots/barplot_taxonomy/Barplots_taxonomy.pdf",
            "results/plots/barplot_qc/Barplot_QC.pdf",
            "results/plots/pca/ACP_results.pdf"
        ]

    # --- DESeq2 + Physico ---
    if config["run_deseq2"]:
        targets += [
            "results/deseq2/deseq2_digesteur_vs_TD1.parquet",
            "results/deseq2/deseq2_date_timeline.parquet",
            "results/deseq2/deseq2_digesteur_combos.parquet"
        ]
    if config["run_physico"]:
        targets.append("results/plots/physico/Physico_plots.pdf")

    return targets

rule all:
    input: get_targets()

# ==========================================================================
# UTILS — Taxonomie NCBI
# ==========================================================================
rule download_taxonomy:
    output:
        zip = temp("data/taxonomy/new_taxdump.zip"),
        csv = "data/taxonomy/mapping_taxons.csv"
    params:
        url      = config["taxonomy"]["zip_url"],
        dmp_name = config["taxonomy"]["dmp_name"]
    conda:   "envs/python_env.yaml"
    shell:
        "wget -q {params.url} -O {output.zip} && "
        "python scripts/utils/download_taxonomy.py"

# ==========================================================================
# UTILS — QC
# ==========================================================================
rule qc_report:
    input:
        contigs = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
        reads   = expand("results/reads/2.Final_Results/{sample}_annotated.tsv",   sample=SAMPLES)
    output:
        tsv = "results/qc/Rapport_QC_Final.tsv",
        pdf = "results/plots/barplot_qc/Barplot_QC.pdf",
        parquet = "Data/Parquet/qc_report.parquet"
    params:
        ordre = config["plots"]["qc"]["ordre"]
    conda:  "envs/r_env.yaml"
    script: "scripts/utils/qc_report.R"

# ==========================================================================
# READS — 2 étapes
# ==========================================================================
rule reads_filter:
    input:
        mapping    = "data/taxonomy/mapping_taxons.csv",
        reads      = lambda w: READS_FILES[w.sample]
    output:
        counts = "results/reads/1.Filtered/{sample}_filtered.tsv"
    params:
        top_n = config["reads"]["top_n"],
        noise = config["reads"]["noise"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/reads/01_filter_reads.py"

rule reads_add_taxaname:
    input:
        data     = "results/reads/1.Filtered/{sample}_filtered.tsv",
        taxonomy = "data/taxonomy/mapping_taxons.csv"
    output:
        taxaname = "results/reads/2.Final_Results/{sample}_annotated.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/reads/02_add_taxaname.py"

# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================
rule contigs_filter_length:
    input:
        raw_data = lambda w: CONTIGS_FILES[w.sample]
    output:
        counts = "results/contigs/1.Raw_Counting/{sample}_counts.tsv"
    params:
        length_threshold = config["contigs"]["length_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/contigs/01_filter_length.py"

rule contigs_filter_abundance:
    input:
        all_samples    = expand("results/contigs/1.Raw_Counting/{sample}_counts.tsv", sample=SAMPLES),
        current_sample = "results/contigs/1.Raw_Counting/{sample}_counts.tsv"
    output:
        counts = "results/contigs/2.Filtered/{sample}_filtered.tsv"
    params:
        abundance_threshold = config["contigs"]["abundance_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/contigs/02_filter_abundance.py"

rule contigs_rpkm:
    input:
        data = "results/contigs/2.Filtered/{sample}_filtered.tsv"
    output:
        rpkm = "results/contigs/3.RPKM/{sample}_rpkm.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/contigs/03_rpkm.py"

rule contigs_rpkm_filter:
    input:
        data           = expand("results/contigs/3.RPKM/{sample}_rpkm.tsv", sample=SAMPLES),
        current_sample = "results/contigs/3.RPKM/{sample}_rpkm.tsv"
    output:
        rpkm = "results/contigs/4.RPKM_Filtered/{sample}_rpkm_filtered.tsv"
    params:
        rpkm_threshold = config["contigs"]["rpkm_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/contigs/04_rpkm_filter.py"

rule contigs_intersection:
    input:
        abund  = "results/contigs/2.Filtered/{sample}_filtered.tsv",
        rpkm   = "results/contigs/4.RPKM_Filtered/{sample}_rpkm_filtered.tsv",
        source = "results/contigs/3.RPKM/{sample}_rpkm.tsv"
    output:
        out = "results/contigs/5.Intersection/{sample}_intersec.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/contigs/05_intersection.py"

rule contigs_add_taxaname:
    input:
        data     = "results/contigs/5.Intersection/{sample}_intersec.tsv",
        taxonomy = "data/taxonomy/mapping_taxons.csv"
    output:
        taxaname = "results/contigs/6.Final_Results/{sample}_annotated.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/contigs/06_add_taxaname.py"

# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================
rule kegg_extraction:
    input:
        raw_data = lambda w: KEGG_FILES[w.sample]
    output:
        counts = "results/kegg/1.Raw_Counting/{sample}_counts.tsv"
    params:
        length_threshold = config["kegg"]["length_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/kegg/01_Kegg_extraction.py"

rule kegg_intersec_count:
    input:
        all_samples    = expand("results/kegg/1.Raw_Counting/{sample}_counts.tsv", sample=SAMPLES),
        current_sample = "results/kegg/1.Raw_Counting/{sample}_counts.tsv"
    output:
        counts = "results/kegg/2.Filtered/{sample}_filtered.tsv"
    params:
        abundance_threshold = config["kegg"]["abundance_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/kegg/02_Kegg_count_intersection.py"

rule kegg_prepared_deseq2:
    input:
        data = "results/kegg/2.Filtered/{sample}_filtered.tsv"
    output:
        rpkm = "results/kegg/3.RPKM/{sample}_rpkm.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/kegg/03_Kegg_prepared_for_DESeq2.py"

rule kegg_standardization:
    input:
        data           = expand("results/kegg/3.RPKM/{sample}_rpkm.tsv", sample=SAMPLES),
        current_sample = "results/kegg/3.RPKM/{sample}_rpkm.tsv"
    output:
        rpkm = "results/kegg/4.RPKM_Filtered/{sample}_rpkm_filtered.tsv"
    params:
        rpkm_threshold = config["kegg"]["rpkm_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "scripts/kegg/021_Kegg_standardization.py"

rule kegg_standardization_agregation:
    input:
        abund  = "results/kegg/2.Filtered/{sample}_filtered.tsv",
        rpkm   = "results/kegg/4.RPKM_Filtered/{sample}_rpkm_filtered.tsv",
        source = "results/kegg/3.RPKM/{sample}_rpkm.tsv"
    output:
        out = "results/kegg/5.Intersection/{sample}_intersec.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/kegg/022_Kegg_standardization_agregation.py"

rule kegg_merge_pathway_levels:
    input:
        data     = "results/kegg/5.Intersection/{sample}_intersec.tsv",
        taxonomy = "data/taxonomy/mapping_taxons.csv"
    output:
        taxaname = "results/kegg/6.Final_Results/{sample}_annotated.tsv"
    conda:   "envs/python_env.yaml"
    script:  "scripts/kegg/023_Kegg_merge_pathway_levels.py"

# ==========================================================================
# PLOTS R
# ==========================================================================

PLOT_SOURCES = {
    "contigs": expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
    "reads":   expand("results/reads/2.Final_Results/{sample}_annotated.tsv",   sample=SAMPLES)
}

PLOT_PARAMS = {
    "contigs": {"value_col": "RPKM",  "label": "RPKM"},
    "reads":   {"value_col": "Reads", "label": "Read counts"}
}

rule plot_heatmap:
    input:
        rds = "Data/RDS/phyloseq.rds",
        tsv = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES)
    output:
        pdf     = "results/plots/heatmap/Heatmaps.pdf",
        parquet = "Data/Parquet/heatmap_data.parquet"
    params:
        top_n        = config["plots"]["top_n"],
        color_opt    = config["plots"]["color_opt"],
        clust_method = config["plots"]["clust_method"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/heatmap.R"

rule plot_barplot:
    input:
        tsv_dir  = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
        metadata = "config/metadata.xlsx"
    output:
        pdf     = "results/plots/barplot/Barplots.pdf",
        parquet = "Data/Parquet/genus_abundance.parquet"
    params:
        top_n      = config["plots"]["top_n"],
        taxon_rank = config["plots"]["target_rank"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/barplot.R"

rule plot_barplot_taxonomy:
    input:
        tsv = lambda w: PLOT_SOURCES[w.source]
    output:
        pdf     = "results/plots/barplot_taxonomy/{source}/Barplots.pdf",
        parquet = "Data/Parquet/barplot_taxonomy_{source}.parquet"
    params:
        target_rank = config["plots"]["target_rank"],
        top_n       = config["plots"]["top_n"],
        value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"],
        label       = lambda w: PLOT_PARAMS[w.source]["label"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/barplot_taxonomy.R"

#rule plot_barplot_taxonomy:
#    input:
#        intersec = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
#        contigs  = expand("results/contigs/3.RPKM/{sample}_rpkm.tsv", sample=SAMPLES)
#    output:
#        pdf     = "results/plots/barplot_taxonomy/Barplots_taxonomy.pdf",
#        parquet = "Data/Parquet/barplot_taxonomy.parquet"
#    params:
#        target_rank = config["plots"]["target_rank"],
#        top_n       = config["plots"]["top_n"]
#    conda:  "envs/r_env.yaml"
#    script: "scripts/plots/barplot_taxonomy.R"

rule plot_pca:
    input:
        parquet  = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
        metadata = "config/metadata.xlsx",
        physico  = "config/physico.xlsx"
    output:
        pdf     = "results/plots/pca/ACP_results.pdf",
        parquet = "Data/Parquet/pca_results.parquet",
        csv     = "results/plots/pca/contributions.csv"
    params:
        dim_x        = config["plots"]["pca"]["dim_x"],
        dim_y        = config["plots"]["pca"]["dim_y"],
        physico_cols = config["plots"]["pca"]["physico_cols"],
        n_top        = config["plots"]["pca"]["n_top"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/pca.R"

rule plot_physico:
    input:
        excel = "config/physico.xlsx"
    output:
        pdf     = "results/plots/physico/Physico_plots.pdf",
        parquet = "Data/Parquet/physico_data.parquet"
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/physico_params.R"

# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================
rule build_phyloseq:
    input:
        csv_dir  = expand("results/reads/2.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
        metadata = "config/metadata.xlsx"
    output:
        rds     = "Data/RDS/phyloseq.rds",
        parquet = "Data/Parquet/phyloseq_long.parquet"
    params:
        top_n = config["plots"]["top_n"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/build_phyloseq.R"

rule run_deseq2:
    input:
        tsv_dir  = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
        metadata = "config/metadata.xlsx"
    output:
        parquet_ref    = "results/deseq2/deseq2_digesteur_vs_TD1.parquet",
        parquet_date   = "results/deseq2/deseq2_date_timeline.parquet",
        parquet_combos = "results/deseq2/deseq2_digesteur_combos.parquet",
        rds_ref        = "Data/RDS/diagdds_digesteur_ref.rds",
        rds_date       = "Data/RDS/diagdds_date.rds",
        rds_combos     = "Data/RDS/diagdds_digesteur_combos.rds"
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/deseq2.R"