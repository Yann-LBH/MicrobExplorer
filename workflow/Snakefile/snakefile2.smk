import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# ==========================================================================
# Échantillons
# ==========================================================================

data           = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES        = data["sample"].tolist()
READS_FILES    = dict(zip(data["sample"], data["reads_file"]))
CONTIGS_FILES  = dict(zip(data["sample"], data["contigs_file"]))
KEGG_FILES     = dict(zip(data["sample"], data["kegg_file"]))

# ==========================================================================
# Cibles finales
# ==========================================================================
def get_targets():
    targets = []

    # --- Taxonomy ---
    if config["run_taxonomy"]:
        targets +=(
            config["path"]["taxonomy_ncbi"]["local_path"], 
            config["path"]["pathway_bakta"]["local_path"]
        ) #####

    # --- Reads ---
    if config["run_reads"]:
        targets.extend(expand(
            "results/treatment/reads/3.Annotated/annotated_{sample}_reads.tsv",
            sample=SAMPLES
        ))

    # --- Contigs ---
    if config["run_contigs"]:
        targets.extend(expand(
            "results/treatment/contigs/6.Annotated/annotated_{sample}_contigs.tsv",
            sample=SAMPLES
        ))

    if config["run_kegg"]:
        targets.extend(expand(
            "results/treatment/kegg/023.Annotated/annotated_{sample}_kegg.tsv",
            sample=SAMPLES
        ))

    # --- QC ---
    # if config["run_qc"]:
    #     targets.append("results/qc/Rapport_QC_Final.pdf")

    # # --- Plots ---
    # if config["run_plots"]:
    #     targets += [
    #         "results/plots/heatmap/Heatmaps.pdf",
    #         expand("results/plots/barplot_taxonomy/{source}/Barplots.pdf",
    #             source=["contigs", "reads"]),
    #         "results/plots/barplot_taxonomy/Barplots_taxonomy.pdf",
    #         "results/plots/barplot_qc/Barplot_QC.pdf",
    #         "results/plots/pca/ACP_results.pdf"
    #     ]

    # # --- DESeq2 + Physico ---
    # if config["run_deseq2"]:
    #     targets += [
    #         "results/deseq2/deseq2_digesteur_vs_TD1.parquet",
    #         "results/deseq2/deseq2_date_timeline.parquet",
    #         "results/deseq2/deseq2_digesteur_combos.parquet"
    #     ]
    # if config["run_physico"]:
    #     targets.append("results/plots/physico/Physico_plots.pdf")

    return targets

rule all:
    input: get_targets()

# ==========================================================================
# UTILS — Taxonomy and pathways
# ==========================================================================
rule download_taxonomy:
    output:
        zip      = temp("data/taxonomy/new_taxdump.zip"),
        dmp_name = config["path"]["taxonomy_ncbi"]["dmp_name"],
        taxonomy = config["path"]["taxonomy_ncbi"]["local_path"]
    params:
        url      = config["path"]["taxonomy_ncbi"]["zip_url"]
    conda:   "envs/python_env.yaml"
    shell:
        "wget -q {params.url} -O {output.zip} && "
        "python scripts/utils/download_NCBI_taxatable.py"

rule download_pathway:
    output:
        tsv = config["path"]["pathway_bakta"]["local_path"]
    params:
        url = config["path"]["pathway_bakta"]["url"]
    conda:   "envs/python_env.yaml"
    shell:
        "wget -q {params.url} -O {output.tsv} && "
        "python scripts/utils/utils_pathway_levels_extraction.py"

# ==========================================================================
# UTILS — QC
# ==========================================================================
rule qc_reads:
    input:
        reads   = expand("results/reads/2.Final_Results/{sample}_annotated.tsv",   sample=SAMPLES)
    output:
        tsv = "report/qc_reads/Rapport_QC_Reads.tsv",
        pdf = "report/qc_reads/Plots_QC_Reads.pdf",
        parquet = "report/qc_reads/Parquet_QC_Reads.parquet"
    conda:  "envs/r_env.yaml"
    script: "../scripts/utils/utils_qc_Reads.R"

rule qc_contigs:
    input:
        contigs = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES)
    output:
        tsv = "report/qc_contigs/Rapport_QC_Contigs.tsv",
        pdf = "report/qc_contigs/Plots_QC_Contigs.pdf",
        parquet = "report/qc_contigs/Parquet_QC_Contigs.parquet"
    conda:  "envs/r_env.yaml"
    script: "../scripts/utils/utils_qc_Contigs.R"

#rule qc_kegg:
#    input:
#        kegg = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES)
#    output:
#        tsv = "report/qc_kegg/Rapport_QC_Kegg.tsv",
#        pdf = "report/qc_kegg/Plots_QC_Kegg.pdf",
#        parquet = "report/qc_kegg/Parquet_QC_Kegg.parquet"
#    conda:  "envs/r_env.yaml"
#    script: "../scripts/utils/utils_qc_Kegg.R"


# ==========================================================================
# READS — 3 étapes
# ==========================================================================
rule reads_countig:
    input:
        raw_data = lambda w: READS_FILES[w.sample]
    output:
        counted = "results/treatment/reads/1.Counted/counted_{sample}_reads.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/01_Reads_counting_raw.py"

rule reads_filter:
    input:
        data     = "results/treatment/reads/1.Counted/counted_{sample}_reads.tsv"
    output:
        filtered = "results/treatment/reads/2.Filtered/filtered_{sample}_reads.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/02_Reads_filter.py"

rule reads_add_taxaname:
    input:
        data     = "results/treatment/reads/2.Filtered/filtered_{sample}_reads.tsv",
        taxonomy = "data/taxonomy_ncbi/taxaname.csv"
    output:
        taxaname = "results/treatment/reads/3.Annotated/annotated_{sample}_reads.tsv"
    params:
        top_n = config["reads"]["top_n"],
        noise = config["reads"]["noise"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/03_Reads_add_taxaname.py"

# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================
rule contigs_counting:
    input:
        raw_data = lambda w: CONTIGS_FILES[w.sample]
    output:
        counted = "results/treatment/contigs/1.Counted/counted_{sample}_contigs.tsv"
    params:
        length_threshold = config["contigs"]["length_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/01_Contigs_counting.py"

rule contigs_filter:
    input:
        data    = "results/treatment/contigs/1.Counted/counted_{sample}_contigs.tsv"
        #current_sample = "results/contigs/1.Counted/counted_{sample}_contigs.tsv"
    output:
        filtered = "results/treatment/contigs/2.Filtered/filtered_{sample}_contigs.tsv"
    params:
        abundance_threshold = config["contigs"]["abundance_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/02_Contigs_filter.py"

rule contigs_rpkm:
    input:
        data = "results/treatment/contigs/2.Filtered/filtered_{sample}_contigs.tsv"
    output:
        rpkm = "results/treatment/contigs/3.RPKM/rpkm_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/03_Contigs_RPKM.py"

rule contigs_rpkm_filter:
    input:
        data = "results/treatment/contigs/3.RPKM/rpkm_{sample}_contigs.tsv"
    output:
        rpkm_filtered = "results/treatment/contigs/4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    params:
        rpkm_threshold = config["contigs"]["rpkm_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/04_Contigs_RPKM_filter.py"

rule contigs_union:
    input:
        abundance  = "results/treatment/contigs/2.Filtered/filtered_{sample}_contigs.tsv",
        rpkm_filtered   = "results/treatment/contigs/4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    output:
        union = "results/treatment/contigs/5.Union/union_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/05_Contigs_union_filtered.py"

rule contigs_add_taxaname:
    input:
        data     = "results/treatment/contigs/5.Union/union_{sample}_contigs.tsv",
        taxonomy = config["path"]["taxonomy_megahit"]
    output:
        taxaname = "results/treatment/contigs/6.Annotated/annotated_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/06_Contigs_add_taxaname.py"

# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================
rule kegg_extraction:
    input:
        raw_data = lambda w: KEGG_FILES[w.sample]
    output:
        extracted = "results/treatment/kegg/1.Extracted/extracted_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/01_Kegg_extraction.py"

rule kegg_intersec_count:
    input:
        data    = "results/treatment/kegg/1.Extracted/extracted_{sample}_kegg.tsv",
        counted = "results/treatment/contigs/1.Counted/counted_{sample}_contigs.tsv" #use contigs counts
    output:
        intersec = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/02_Kegg_count_intersection.py"

rule kegg_prepared_deseq2:
    input:
        data = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        deseq2 = "results/treatment/kegg/3.Deseq2/deseq2_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/03_Kegg_prepared_for_DESeq2.py"

rule kegg_standardization:
    input:
        data           = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv",
        current_sample = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        stand = "results/treatment/kegg/021.Standardized/standardized_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/021_Kegg_standardization.py"

rule kegg_standardization_agregation:
    input:
        data  = "results/treatment/kegg/021.Standardized/standardized_{sample}_kegg.tsv"
    output:
        agreg = "results/treatment/kegg/022.Aggregated/stand_aggreg_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/022_Kegg_standardization_aggregation.py"

rule kegg_merge_pathway_levels:
    input:
        data     = "results/treatment/kegg/022.Aggregated/stand_aggreg_{sample}_kegg.tsv",
        pathway = config["path"]["pathway_bakta"]["local_path"]
    output:
        taxaname = "results/treatment/kegg/023.Annotated/annotated_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/023_Kegg_merge_pathway_levels.py"

# ==========================================================================
# PLOTS R
# ==========================================================================

# PLOT_SOURCES = {
#     "reads":   expand("results/reads/2.Final_Results/{sample}_annotated.tsv",   sample=SAMPLES),
#     "contigs": expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
#     "kegg": expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES)
# }

# PLOT_PARAMS = {
#     "reads":   {"value_col": "Reads", "label": "Read counts"},
#     "contigs": {"value_col": "RPKM",  "label": "RPKM"},
#     "kegg": {"value_col": "RPKM",  "label": "RPKM"}
# }

# rule plot_heatmap:
#     input:
#         rds = "Data/RDS/phyloseq.rds",
#         data = expand("results/contigs/6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES)
#     output:
#         pdf     = "results/plots/heatmap/Heatmaps.pdf",
#         parquet = "Data/Parquet/heatmap_data.parquet"
#     params:
#         shared       = config["plots"]["shared"],
#         top_n        = config["plots"]["heatmap"]["top_n"],
#         clust_method = config["plots"]["heatmap"]["clust_method"]
#     conda:  "envs/r_env.yaml"
#     script: "scripts/plots/heatmap.R"

# rule plot_barplot:
#     input:
#         tsv_dir  = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
#         metadata = config["path"]["metadata"]
#     output:
#         pdf     = "results/plots/barplot/Barplots.pdf",
#         parquet = "Data/Parquet/genus_abundance.parquet"
#     params:
#         shared     = config["plots"]["shared"],
#         top_n      = config["plots"]["barplot"]["top_n"],
#         taxon_rank = config["plots"]["barplot"]["taxon_rank"]
#     conda:  "envs/r_env.yaml"
#     script: "scripts/plots/barplot.R"

# rule plot_barplot_taxonomy:
#     input:
#         tsv = lambda w: PLOT_SOURCES[w.source]
#     output:
#         pdf     = "results/plots/barplot_taxonomy/{source}/Barplots.pdf",
#         parquet = "Data/Parquet/barplot_taxonomy_{source}.parquet"
#     params:
#         shared      = config["plots"]["shared"],
#         top_n       = config["plots"]["barplot_taxonomy"]["top_n"],
#         taxon_rank  = config["plots"]["barplot_taxonomy"]["taxon_rank"],
#         value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"],
#         label       = lambda w: PLOT_PARAMS[w.source]["label"]
#     conda:  "envs/r_env.yaml"
#     script: "scripts/plots/barplot_taxonomy.R"

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

# rule plot_pca:
#     input:
#         parquet  = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
#         metadata = config["path"]["metadata"],
#         physico  = config["path"]["physico_params"]
#     output:
#         pdf     = "results/plots/pca/ACP_results.pdf",
#         parquet = "Data/Parquet/pca_results.parquet",
#         csv     = "results/plots/pca/contributions.csv"
#     params:
#         shared       = config["plots"]["shared"],
#         top_n        = config["plots"]["pca"]["top_n"],
#         point_size   = config["plots"]["pca"]["point_size"],
#         dim_x        = config["plots"]["pca"]["dim_x"],
#         dim_y        = config["plots"]["pca"]["dim_y"],
#         physico_cols = config["plots"]["pca"]["physico_cols"],
#     conda:  "envs/r_env.yaml"
#     script: "../scripts/plots/pca.R"

# rule plot_physico:
#     input:
#         physico  = config["path"]["physico_params"]
#     output:
#         pdf     = "results/plots/physico/physico_plots.pdf",
#         parquet = "results/parquet/physico_data.parquet"
#     conda:  "envs/r_env.yaml"
#     script: "../scripts/plots/physico_params.R"

# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================
# rule run_phyloseq:
#     input:
#         files_dir  = expand("results/reads/3.Annotated/annotated_{sample}_reads.tsv", sample=SAMPLES),
#         metadata = config["path"]["metadata"]
#     output:
#         rds_reads     = config["save"]["rds_reads"] + "phyloseq_reads.rds",
#         rds_contigs   = config["save"]["rds_contigs"] + "phyloseq_contigs.rds",
#         rds_kegg      = config["save"]["rds_kegg"] + "phyloseq_kegg.rds"
#     conda:  "envs/r_env.yaml"
#     script: "../scripts/plots/build_phyloseq.R"

# rule run_deseq2:
#     input:
#         tsv_dir  = expand("results/contigs/6.Final_Results/{sample}_annotated.tsv", sample=SAMPLES),
#         metadata = config["path"]["metadata"]
#     output:
#         parquet_ref    = "results/deseq2/deseq2_digesteur_vs_TD1.parquet",
#         parquet_date   = "results/deseq2/deseq2_date_timeline.parquet",
#         parquet_combos = "results/deseq2/deseq2_digesteur_combos.parquet",
#         rds_ref        = "Data/RDS/diagdds_digesteur_ref.rds",
#         rds_date       = "Data/RDS/diagdds_date.rds",
#         rds_combos     = "Data/RDS/diagdds_digesteur_combos.rds"
#     conda:  "envs/r_env.yaml"
#     script: "../scripts/plots/deseq2.R"