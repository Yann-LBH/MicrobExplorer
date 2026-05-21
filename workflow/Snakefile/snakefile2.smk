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

    # --- Plots ---
    if config["run_plots"]:
        targets += [
            "results/plots/heatmap/Heatmap_contigs.pdf",
            # Le '*' extrait les éléments de la liste expand() directement ici
            *expand("results/plots/stackedbarplot/{source}/Stackedbarplot_deseq2_{source}.pdf",
                source=["contigs", "kegg"]),
            *expand("results/plots/stackedbarplot/{source}/Stackedbarplot_organisms_abundance_{source}.pdf",
                source=["reads", "contigs"]),
            "results/plots/stackedbarplot/{source}/Stackedbarplot_pathway_abundance_{source}.pdf",
            *expand("results/plots/stackedbarplot/{source}/Stackedbarplot_taxonomy_{source}.pdf",
                source=["reads", "contigs"]),
            *expand("results/plots/pca/{source}/PCA_{source}.pdf",
                source=["contigs", "kegg"]),
            *expand("results/plots/volcano/{source}/Volcano_deseq2_{source}.pdf",
                source=["contigs", "kegg"])
        ]

    # --- DESeq2 + Physico ---
    if config["run_deseq2"]:
        datatypes = ["contigs", "kegg"]
        contrasts = ["ref", "date", "combo"]

        targets += [
            config["save"][f"parquet_{d}"] + f"deseq2_{c}_{d}.parquet"
            for d in datatypes
            for c in contrasts
        ]
        targets += [
            config["save"][f"rds_{d}"] + f"deseq2_{c}_{d}.rds"
            for d in datatypes
            for c in contrasts
        ]

    #if config["run_physico"]:
    #    targets.append("results/plots/physico/Physico_plots.pdf")
    
    # --- QC ---
    # if config["run_qc"]:
    #     targets.append("results/qc/Rapport_QC_Final.pdf")
  
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
        reads       = "results/reads/2.Final_Results/{sample}_annotated.tsv",
    output:
        tsv         = "report/qc_reads/Rapport_QC_Reads.tsv",
        pdf         = "report/qc_reads/Plots_QC_Reads.pdf",
        parquet     = "report/qc_reads/Parquet_QC_Reads.parquet"
    conda:  "envs/r_env.yaml"
    script: "../scripts/utils/utils_qc_Reads.R"

rule qc_contigs:
    input:
        contigs     = "results/contigs/6.Final_Results/{sample}_annotated.tsv",
    output:
        tsv         = "report/qc_contigs/Rapport_QC_Contigs.tsv",
        pdf         = "report/qc_contigs/Plots_QC_Contigs.pdf",
        parquet     = "report/qc_contigs/Parquet_QC_Contigs.parquet"
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
        raw_data    = lambda w: READS_FILES[w.sample]
    output:
        counted     = "results/treatment/reads/1.Counted/counted_{sample}_reads.tsv"
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
        data        = "results/treatment/reads/2.Filtered/filtered_{sample}_reads.tsv",
        taxonomy    = "data/taxonomy_ncbi/taxaname.csv"
    output:
        taxaname    = "results/treatment/reads/3.Annotated/annotated_{sample}_reads.tsv"
    params:
        top_n       = config["reads"]["top_n"],
        noise       = config["reads"]["noise"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/03_Reads_add_taxaname.py"

# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================
rule contigs_counting:
    input:
        raw_data            = lambda w: CONTIGS_FILES[w.sample]
    output:
        counted             = "results/treatment/contigs/1.Counted/counted_{sample}_contigs.tsv"
    params:
        length_threshold    = config["contigs"]["length_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/01_Contigs_counting.py"

rule contigs_filter:
    input:
        data                = "results/treatment/contigs/1.Counted/counted_{sample}_contigs.tsv"
        #current_sample = "results/contigs/1.Counted/counted_{sample}_contigs.tsv"
    output:
        filtered            = "results/treatment/contigs/2.Filtered/filtered_{sample}_contigs.tsv"
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
        data            = "results/treatment/contigs/3.RPKM/rpkm_{sample}_contigs.tsv"
    output:
        rpkm_filtered   = "results/treatment/contigs/4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    params:
        rpkm_threshold  = config["contigs"]["rpkm_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/04_Contigs_RPKM_filter.py"

rule contigs_union:
    input:
        abundance       = "results/treatment/contigs/2.Filtered/filtered_{sample}_contigs.tsv",
        rpkm_filtered   = "results/treatment/contigs/4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    output:
        union           = "results/treatment/contigs/5.Union/union_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/05_Contigs_union_filtered.py"

rule contigs_add_taxaname:
    input:
        data        = "results/treatment/contigs/5.Union/union_{sample}_contigs.tsv",
        taxonomy    = config["path"]["taxonomy_megahit"]
    output:
        taxaname    = "results/treatment/contigs/6.Annotated/annotated_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/06_Contigs_add_taxaname.py"

# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================
rule kegg_extraction:
    input:
        raw_data    = lambda w: KEGG_FILES[w.sample]
    output:
        extracted   = "results/treatment/kegg/1.Extracted/extracted_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/01_Kegg_extraction.py"

rule kegg_intersec_count:
    input:
        data        = "results/treatment/kegg/1.Extracted/extracted_{sample}_kegg.tsv",
        counted     = "results/treatment/contigs/1.Counted/counted_{sample}_contigs.tsv" #use contigs counts
    output:
        intersec    = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/02_Kegg_count_intersection.py"

rule kegg_prepared_deseq2:
    input:
        data    = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        deseq2  = "results/treatment/kegg/3.Deseq2/deseq2_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/03_Kegg_prepared_for_DESeq2.py"

rule kegg_standardization:
    input:
        data            = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv",
        current_sample  = "results/treatment/kegg/2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        stand           = "results/treatment/kegg/021.Standardized/standardized_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/021_Kegg_standardization.py"

rule kegg_standardization_agregation:
    input:
        data    = "results/treatment/kegg/021.Standardized/standardized_{sample}_kegg.tsv"
    output:
        agreg   = "results/treatment/kegg/022.Aggregated/stand_aggreg_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/022_Kegg_standardization_aggregation.py"

rule kegg_merge_pathway_levels:
    input:
        data        = "results/treatment/kegg/022.Aggregated/stand_aggreg_{sample}_kegg.tsv",
        pathway     = config["path"]["pathway_bakta"]["local_path"]
    output:
        taxaname    = "results/treatment/kegg/023.Annotated/annotated_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/023_Kegg_merge_pathway_levels.py"

# ==========================================================================
# PLOTS R
# ==========================================================================

PLOT_SOURCES = {
    "reads":    expand("results/treatment/reads/3.Annotated/annotated_{sample}_reads.tsv",   sample=SAMPLES),
    "contigs":  expand("results/treatment/contigs/6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES),
    "kegg":     expand("results/treatment/kegg/023.Annotated/annotated_{sample}_kegg.tsv", sample=SAMPLES)
}

PLOT_PARAMS = {
    "reads":    {"value_col": "Reads", "label": "Read counts"},
    "contigs":  {"value_col": "RPKM",  "label": "RPKM"},
    "kegg":     {"value_col": "RPKM",  "label": "RPKM"}
}

rule plot_stackedbarplot_DESeq2:
    input:
        data        = lambda w: PLOT_SOURCES[w.source],
        metadata    = config["path"]["metadata"]
    output:
        pdf         = "results/plots/stackedbarplot/{source}/Stackedbarplot_deseq2_{source}.pdf",
        parquet     = "results/parquet/stackedbarplot_deseq2_{source}.parquet"
    params:
        shared      = config["plots"]["shared"],
        top_n       = config["plots"]["barplot"]["top_n"],
        taxon_rank  = config["plots"]["barplot"]["taxon_rank"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/Stackedbarplot_from_DESeq2.R"

rule plot_stackedbarplot_organisms:
    input:
        data        = lambda w: PLOT_SOURCES[w.source],
        metadata    = config["path"]["metadata"]
    output:
        pdf         = "results/plots/stackedbarplot/{source}/Stackedbarplot_organisms_abundance_{source}.pdf",
        parquet     = "results/parquet/stackedbarplot_organisms_abundance_{source}.parquet"
    params:
        shared      = config["plots"]["shared"],
        top_n       = config["plots"]["stackedbarplot"]["top_n"],
        taxon_rank  = config["plots"]["stackedbarplot"]["taxon_rank"],
        value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"],
        label       = lambda w: PLOT_PARAMS[w.source]["label"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/Stackedbarplot_organisms_abundance.R"

rule plot_stackedbarplot_pathway:
    input:
        data        = lambda w: PLOT_SOURCES["kegg"],
        metadata    = config["path"]["metadata"]
    output:
        pdf         = "results/plots/stackedbarplot/kegg/Stackedbarplot_pathway_abundance_kegg.pdf",
        parquet     = "results/parquet/stackedbarplot_pathway_abundance_kegg.parquet"
    params:
        shared      = config["plots"]["shared"],
        top_n       = config["plots"]["stackedbarplot"]["top_n"],
        taxon_rank  = config["plots"]["stackedbarplot"]["taxon_rank"],
        value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"],
        label       = lambda w: PLOT_PARAMS[w.source]["label"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/Stackedbarplot_pathway_abundance.R"

rule plot_stackedbarplot_taxonomy:
    input:
        data        = lambda w: PLOT_SOURCES[w.source],
        metadata    = config["path"]["metadata"]
        #intersec    = expand("results/treatment/contigs/6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES),
        #contigs     = expand("results/treatment/contigs/3.RPKM/rpkm_{sample}_contigs.tsv", sample=SAMPLES)
    output:
        pdf         = "results/plots/stackedbarplot/{source}/Stackedbarplots_taxonomy_{source}.pdf",
        parquet     = "results/parquet/stackedbarplot_taxonomy_{source}.parquet"
    params:
        taxon_rank  = config["plots"]["stackedbarplot"]["taxon_rank"],
        top_n       = config["plots"]["stackedbarplot"]["top_n"]
    conda:  "envs/r_env.yaml"
    script: "scripts/plots/Stackedbarplot_taxonomy.R"

rule plot_heatmap: # CHECK
    input:
        rds             = config["save"]["rds_contigs"] + "phyloseq_contigs.rds",
        data            = PLOT_SOURCES["contigs"],
        metadata        = config["path"]["metadata"]
    output:
        pdf             = "results/plots/heatmap/Heatmap_contigs.pdf",
        parquet         = config["save"]["parquet_contigs"] + "heatmap_contigs.parquet",
    params:
        shared          = config["plots"]["shared"],
        top_n           = config["plots"]["heatmap"]["top_n"],
        clust_method    = config["plots"]["heatmap"]["clust_method"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Heatmap.R"

rule plot_pca:
    input:
        data            = lambda w: PLOT_SOURCES[w.source],
        metadata        = config["path"]["metadata"],
        physico         = config["path"]["physico_params"]
    output:
        pdf             = "results/plots/pca/{source}/PCA_{source}.pdf",
        parquet         = "Data/Parquet/pca_results.parquet",
        csv             = "results/plots/pca/contributions.csv"
    params:
        shared          = config["plots"]["shared"],
        top_n           = config["plots"]["pca"]["top_n"],
        point_size      = config["plots"]["pca"]["point_size"],
        dim_x           = config["plots"]["pca"]["dim_x"],
        dim_y           = config["plots"]["pca"]["dim_y"],
        physico_cols    = config["plots"]["pca"]["physico_cols"],
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/PCA.R"

rule plot_volcano_DESeq2:
    input:
        rds         = config["save"]["rds_contigs"] + "deseq2_contigs.rds",
        data        = lambda w: PLOT_SOURCES[w.source]
    output:
        pdf         = "results/plots/volcano/{source}/Volcano_deseq2_{source}.pdf",
        parquet     = "results/parquet/volcano_from_deseq2.parquet"
    params:
        padj        = config["plots"]["run_volcano"]["pvalue_threshold"],
        lfc         = config["plots"]["run_volcano"]["lfc_treshold"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Volcano_from_DESeq2.R"

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

ANALYSIS_SOURCES = {
    "reads":        expand("results/treatment/reads/3.Annotated/annotated_{sample}_reads.tsv",   sample=SAMPLES),
    "contigs":      expand("results/treatment/contigs/6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES),
    "kegg":         expand("results/treatment/kegg/023.Annotated/annotated_{sample}_kegg.tsv", sample=SAMPLES),
    "kegg_deseq2":  expand("results/treatment/kegg/3.Deseq2/deseq2_{sample}_kegg.tsv", sample=SAMPLES)
}

rule run_phyloseq:
    input:
        reads_dir       = ANALYSIS_SOURCES["reads"],
        contigs_dir     = ANALYSIS_SOURCES["contigs"],
        kegg_dir        = ANALYSIS_SOURCES["kegg"],
        metadata        = config["path"]["metadata"]
    output:
        rds_reads       = config["save"]["rds_reads"] + "phyloseq_reads.rds",
        rds_contigs     = config["save"]["rds_contigs"] + "phyloseq_contigs.rds",
        rds_kegg        = config["save"]["rds_kegg"] + "phyloseq_kegg.rds"
    conda:  "envs/r_env.yaml"
    script: "../scripts/analysis/Phyloseq.R"

rule run_deseq2:
    input:
        contigs_dir         = ANALYSIS_SOURCES["contigs"],
        kegg_dir            = ANALYSIS_SOURCES["kegg_deseq2"],
        metadata            = config["path"]["metadata"]
    output:
        # --- OBJETS COMPLETS (RDS) ---
        rds_ref_r           = config["save"]["rds_reads"] + "deseq2_ref_reads.rds",
        rds_date_r          = config["save"]["rds_reads"] + "deseq2_date_reads.rds",
        rds_combo_r         = config["save"]["rds_reads"] + "deseq2_combo_reads.rds",
        rds_ref_c           = config["save"]["rds_contigs"] + "deseq2_ref_contigs.rds",
        rds_date_c          = config["save"]["rds_contigs"] + "deseq2_date_contigs.rds",
        rds_combo_c         = config["save"]["rds_contigs"] + "deseq2_combo_contigs.rds",
        rds_ref_k           = config["save"]["rds_kegg"] + "deseq2_ref_kegg.rds",
        rds_date_k          = config["save"]["rds_kegg"] + "deseq2_date_kegg.rds",
        rds_combo_k         = config["save"]["rds_kegg"] + "deseq2_combo_kegg.rds",

        # --- TABLES DE RÉSULTATS (PARQUET) ---
        parquet_ref_r       = config["save"]["parquet_reads"] + "deseq2_ref_reads.parquet",
        parquet_date_r      = config["save"]["parquet_reads"] + "deseq2_date_reads.parquet",
        parquet_combo_r     = config["save"]["parquet_reads"] + "deseq2_combo_reads.parquet",
        parquet_ref_c       = config["save"]["parquet_contigs"] + "deseq2_ref_contigs.parquet",
        parquet_date_c      = config["save"]["parquet_contigs"] + "deseq2_date_contigs.parquet",
        parquet_combo_c     = config["save"]["parquet_contigs"] + "deseq2_combo_contigs.parquet",
        parquet_ref_k       = config["save"]["parquet_kegg"] + "deseq2_ref_kegg.parquet",
        parquet_date_k      = config["save"]["parquet_kegg"] + "deseq2_date_kegg.parquet",
        parquet_combo_k     = config["save"]["parquet_kegg"] + "deseq2_combo_kegg.parquet"

    conda:  "envs/r_env.yaml"
    script: "../scripts/analysis/DESeq2.R"