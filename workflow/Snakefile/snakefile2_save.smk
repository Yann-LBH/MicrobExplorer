import pandas as pd
from input_pathlib import input_path

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
            config["input_path"]["taxonomy_ncbi"]["local_input_path"], 
            config["input_path"]["input_pathway_bakta"]["local_input_path"]
        ) #####

    # --- Reads ---
    if config["run_reads"]:
        targets.extend(expand(
            READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv",
            sample=SAMPLES
        ))

    # --- Contigs ---
    if config["run_contigs"]:
        targets.extend(expand(
            CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv",
            sample=SAMPLES
        ))

    if config["run_kegg"]:
        targets.extend(expand(
            KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv",
            sample=SAMPLES
        ))

    MODE_TO_FILENAME = {
        "pathway":   "pathway_abundance",
        "organisms": "organisms_abundance",
        "taxonomy":  "taxonomy",
    }
    
    # --- Plots ---
    if config["run_plots"]:
        targets += [
            # Heatmap
            *expand(config["output_path"]["plots"] + "heatmap/Heatmap_{source}.pdf",
                source=config["datatypes"]["heatmap"]),
            
            # Stackedbarplots
            *expand(config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplot_deseq2_{source}.pdf",
                source=config["datatypes"]["stackedbarplot"]["deseq2"]),
                
            *expand(config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplot_organisms_abundance_{source}.pdf",
                source=config["datatypes"]["stackedbarplot"]["organisms"]),
                
            *expand(config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplot_input_pathway_abundance_{source}.pdf",
                source=config["datatypes"]["stackedbarplot"]["pathway"]),
                
            *expand(config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplot_taxonomy_{source}.pdf",
                source=config["datatypes"]["stackedbarplot"]["taxonomy"]),
            for cle, sources in config["datatypes"]["stackedbarplot"].items():
        
            # Pour pathway, le nom du fichier est un peu différent ("input_pathway_abundance")
            # On ajuste le nom du fichier dynamiquement selon la clé
            if cle == "pathway":
                nom_fichier = "input_pathway_abundance"
            elif cle == "organisms":
                nom_fichier = "organisms_abundance"
            else:
                nom_fichier = cle  # Sera "deseq2" ou "taxonomy"
                
            # On ajoute dynamiquement les cibles générées par expand
            targets += expand(
                config["output_path"]["plots"] + f"stackedbarplot/{{source}}/Stackedbarplot_{nom_fichier}_{{source}}.pdf",
                source=sources
            )
            
            # Volcano
            *expand(config["output_path"]["plots"] + "volcano/{source}/Volcano_deseq2_{source}.pdf",
                source=config["datatypes"]["volcano"])
        ]

    if config["run_phyloseq"]:
        datatypes = ["reads", "contigs", "kegg"]

        targets += [
            config["output_path"]["rds"] + f"{d}/phyloseq_{d}.rds"
            for d in datatypes
        ]

        targets += [
            config["output_path"]["parquet"] + f"{d}/phyloseq_{d}.parquet"
            for d in datatypes
        ]

    # --- DESeq2 ---
    if config["run_deseq2"]:
        datatypes = ["contigs", "kegg"]
        contrasts = ["ref", "date", "combo"]

        targets += [
            config["output_path"]["rds"] + f"{d}/deseq2_{c}_{d}.rds"
            for d in datatypes
            for c in contrasts
        ]

        targets += [
            config["output_path"]["parquet"] + f"{d}/deseq2_{c}_{d}.parquet"
            for d in datatypes
            for c in contrasts
        ]

    #if config["run_physico"]:
    #    targets.append(["output_path"]["plots"] + "physico/Physico_plots.pdf")
    
    # --- QC ---
    # if config["run_qc"]:
    #     targets.append("results/qc/Rapport_QC_Final.pdf")
  
    return targets

rule all:
    input: get_targets()

# ==========================================================================
# UTILS — Taxonomy and input_pathways
# ==========================================================================
rule download_taxonomy:
    output:
        zip      = temp("data/taxonomy/new_taxdump.zip"),
        dmp_name = config["input_path"]["taxonomy_ncbi"]["dmp_name"],
        taxonomy = config["input_path"]["taxonomy_ncbi"]["local_input_path"]
    params:
        url      = config["input_path"]["taxonomy_ncbi"]["zip_url"]
    conda:   "envs/python_env.yaml"
    shell:
        "wget -q {params.url} -O {output.zip} && "
        "python scripts/utils/download_NCBI_taxatable.py"

rule download_input_pathway:
    output:
        tsv = config["input_path"]["input_pathway_bakta"]["local_input_path"]
    params:
        url = config["input_path"]["input_pathway_bakta"]["url"]
    conda:   "envs/python_env.yaml"
    shell:
        "wget -q {params.url} -O {output.tsv} && "
        "python scripts/utils/utils_input_pathway_levels_extraction.py"

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
READS_TREATMENT = config["output_path"]["treatment"] + "reads/"

rule reads_countig:
    input:
        raw_data    = lambda w: READS_FILES[w.sample]
    output:
        counted     = READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/01_Reads_counting_raw.py"

rule reads_filter:
    input:
        data     = READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    output:
        filtered = READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/02_Reads_filter.py"

rule reads_add_taxaname:
    input:
        data        = READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv",
        taxonomy    = config["input_path"]["taxonomy_ncbi"]["local_path"]
    output:
        taxaname    = READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv"
    params:
        top_n       = config["reads"]["top_n"],
        noise       = config["reads"]["noise"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/reads/03_Reads_add_taxaname.py"

# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================
CONTIGS_TREATMENT = config["output_path"]["treatment"] + "contigs/"

rule contigs_counting:
    input:
        raw_data            = lambda w: CONTIGS_FILES[w.sample]
    output:
        counted             = CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    params:
        length_threshold    = config["contigs"]["length_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/01_Contigs_counting.py"

rule contigs_filter:
    input:
        data                = CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
        #current_sample = "results/contigs/1.Counted/counted_{sample}_contigs.tsv"
    output:
        filtered            = CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv"
    params:
        abundance_threshold = config["contigs"]["abundance_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/02_Contigs_filter.py"

rule contigs_rpkm:
    input:
        data = CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv"
    output:
        rpkm = CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/03_Contigs_RPKM.py"

rule contigs_rpkm_filter:
    input:
        data            = CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv"
    output:
        rpkm_filtered   = CONTIGS_TREATMENT + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    params:
        rpkm_threshold  = config["contigs"]["rpkm_threshold"]
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/04_Contigs_RPKM_filter.py"

rule contigs_union:
    input:
        abundance       = CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv",
        rpkm_filtered   = CONTIGS_TREATMENT + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    output:
        union           = CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/05_Contigs_union_filtered.py"

rule contigs_add_taxaname:
    input:
        data        = CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv",
        taxonomy    = config["input_path"]["taxonomy_megahit"]
    output:
        taxaname    = CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/contigs/06_Contigs_add_taxaname.py"

# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================
KEGG_TREATMENT = config["output_path"]["treatment"] + "kegg/"

rule kegg_extraction:
    input:
        raw_data    = lambda w: KEGG_FILES[w.sample]
    output:
        extracted   = KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/01_Kegg_extraction.py"

rule kegg_intersec_count:
    input:
        data        = KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv",
        counted     = CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv" #use contigs counts
    output:
        intersec    = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/02_Kegg_count_intersection.py"

rule kegg_prepared_deseq2:
    input:
        data    = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        deseq2  = KEGG_TREATMENT + "3.Deseq2/deseq2_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/03_Kegg_prepared_for_DESeq2.py"

rule kegg_standardization:
    input:
        data            = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv",
        current_sample  = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        stand           = KEGG_TREATMENT + "021.Standardized/standardized_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/021_Kegg_standardization.py"

rule kegg_standardization_agregation:
    input:
        data    = KEGG_TREATMENT + "021.Standardized/standardized_{sample}_kegg.tsv"
    output:
        agreg   = KEGG_TREATMENT + "022.Aggregated/stand_aggreg_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/022_Kegg_standardization_aggregation.py"

rule kegg_merge_input_pathway_levels:
    input:
        data        = KEGG_TREATMENT + "022.Aggregated/stand_aggreg_{sample}_kegg.tsv",
        input_pathway     = config["input_path"]["input_pathway_bakta"]["local_input_path"]
    output:
        taxaname    = KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv"
    conda:   "envs/python_env.yaml"
    script:  "../scripts/kegg/023_Kegg_merge_input_pathway_levels.py"

# ==========================================================================
# PLOTS R
# ==========================================================================

TREATMENT_SOURCES = {
    "reads":        expand(READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv",   sample=SAMPLES),
    "contigs":      expand(CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES),
    "kegg":         expand(KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv", sample=SAMPLES),
    "kegg_deseq2":  expand(KEGG_TREATMENT + "3.Deseq2/deseq2_{sample}_kegg.tsv", sample=SAMPLES)
}

PLOT_PARAMS = {
    "reads":    {"value_col": "Reads", "label": "Read counts"},
    "contigs":  {"value_col": "RPKM",  "label": "RPKM"},
    "kegg":     {"value_col": "RPKM",  "label": "RPKM"}
}

rule plot_stackedbarplot_DESeq2:
    input:
        data        = lambda w: TREATMENT_SOURCES[w.source],
        metadata    = config["input_path"]["metadata"]
    output:
        pdf         = config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplot_deseq2_{source}.pdf",
        parquet     = config["output_path"]["parquet"] + "{source}/stackedbarplot_deseq2_{source}.parquet"
    params:
        shared      = config["plots"]["shared"],
        top_n       = config["plots"]["barplot"]["top_n"],
        taxon_rank  = config["plots"]["barplot"]["taxon_rank"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Stackedbarplot_from_DESeq2.R"

rule plot_stackedbarplot_organisms:
    input:
        data        = lambda w: TREATMENT_SOURCES[w.source],
        metadata    = config["input_path"]["metadata"]
    output:
        pdf         = config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplot_organisms_abundance_{source}.pdf",
        parquet     = ["output_path"]["parquet"] + "{source}/stackedbarplot_organisms_abundance_{source}.parquet"
    params:
        shared      = config["plots"]["shared"],
        top_n       = config["plots"]["stackedbarplot"]["top_n"],
        taxon_rank  = config["plots"]["stackedbarplot"]["taxon_rank"],
        value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"],
        label       = lambda w: PLOT_PARAMS[w.source]["label"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Stackedbarplot_organisms_abundance.R"

rule plot_stackedbarplot_input_pathway:
    input:
        data        = lambda w: TREATMENT_SOURCES["kegg"],
        metadata    = config["input_path"]["metadata"]
    output:
        pdf         = config["output_path"]["plots"] + "stackedbarplot/kegg/Stackedbarplot_input_pathway_abundance_kegg.pdf",
        parquet     = config["output_path"]["parquet"] + "{source}/stackedbarplot_input_pathway_abundance_kegg.parquet"
    params:
        shared      = config["plots"]["shared"],
        top_n       = config["plots"]["stackedbarplot"]["top_n"],
        taxon_rank  = config["plots"]["stackedbarplot"]["taxon_rank"],
        value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"],
        label       = lambda w: PLOT_PARAMS[w.source]["label"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Stackedbarplot_input_pathway_abundance.R"

rule plot_stackedbarplot_taxonomy:
    input:
        data        = lambda w: TREATMENT_SOURCES[w.source],
        metadata    = config["input_path"]["metadata"]
        #intersec    = expand(CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES),
        #contigs     = expand(CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv", sample=SAMPLES)
    output:
        pdf         = config["output_path"]["plots"] + "stackedbarplot/{source}/Stackedbarplots_taxonomy_{source}.pdf",
        parquet     = config["output_path"]["parquet"] + "{source}/stackedbarplot_taxonomy_{source}.parquet"
    params:
        taxon_rank  = config["plots"]["stackedbarplot"]["taxon_rank"],
        top_n       = config["plots"]["stackedbarplot"]["top_n"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Stackedbarplot_taxonomy.R"

rule plot_heatmap: # CHECK
    input:
        rds             = config["output_path"]["rds_contigs"] + "phyloseq_contigs.rds",
        data            = TREATMENT_SOURCES["contigs"],
        metadata        = config["input_path"]["metadata"]
    output:
        pdf             = config["output_path"]["plots"] + "contigs/heatmap/Heatmap_contigs.pdf",
        parquet         = config["output_path"]["parquet"] + "contigs/heatmap_contigs.parquet",
    params:
        shared          = config["plots"]["shared"],
        top_n           = config["plots"]["heatmap"]["top_n"],
        clust_method    = config["plots"]["heatmap"]["clust_method"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Heatmap.R"

rule plot_pca:
    input:
        data            = lambda w: TREATMENT_SOURCES[w.source],
        metadata        = config["input_path"]["metadata"],
        physico         = config["input_path"]["physico_params"]
    output:
        pdf             = config["output_path"]["plots"] + "pca/{source}/PCA_{source}.pdf",
        parquet         = config["output_path"]["parquet"] + "{source}/pca_{source}.parquet",
        csv             = "results/infos/pca/pca_contributions_{source}.csv"
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
        rds         = config["output_path"]["rds_contigs"] + "deseq2_contigs.rds",
        data        = lambda w: TREATMENT_SOURCES[w.source]
    output:
        pdf         = config["output_path"]["plots"] + "volcano/{source}/Volcano_deseq2_{source}.pdf",
        parquet     = config["output_path"]["parquet"] + "{source}/volcano_from_deseq2_{source}.parquet"
    params:
        padj        = config["plots"]["run_volcano"]["pvalue_threshold"],
        lfc         = config["plots"]["run_volcano"]["lfc_treshold"]
    conda:  "envs/r_env.yaml"
    script: "../scripts/plots/Volcano_from_DESeq2.R"

# rule plot_physico:
#     input:
#         physico  = config["input_path"]["physico_params"]
#     output:
#         pdf     = config["output_path"]["plots"] + "physico/physico_plots.pdf",
#         parquet = config["output_path"]["parquet"] + "{source}/physico_data.parquet"
#     conda:  "envs/r_env.yaml"
#     script: "../scripts/plots/physico_params.R"

# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================

rule run_phyloseq:
    input:
        data        = lambda w: TREATMENT_SOURCES[w.d],
        metadata    = config["input_path"]["metadata"]
    output:
        rds         = config["output_path"]["rds"] + "{d}/phyloseq_{d}.rds",
    conda:  "envs/r_env.yaml"
    script: "../scripts/analysis/Phyloseq.R"

rule run_deseq2:
    input:
        data        = lambda w: TREATMENT_SOURCES["kegg_deseq2" if w.d == "kegg" else w.d],
        metadata    = config["input_path"]["metadata"]
    output:
        rds         = config["output_path"]["rds"] + "{d}/deseq2_{c}_{d}.rds",
        parquet     = config["output_path"]["parquet"] + "{d}/deseq2_{c}_{d}.parquet"

    conda:  "envs/r_env.yaml"
    script: "../scripts/analysis/DESeq2.R"