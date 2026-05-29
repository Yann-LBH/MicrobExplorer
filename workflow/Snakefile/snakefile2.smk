import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# ==========================================================================
#   DICTIONARIES
# ==========================================================================

# Échantillons
data           = pd.read_csv(config["samples_file"], sep="\t")
SAMPLES        = data["sample"].tolist()
READS_FILES    = dict(zip(data["sample"], data["reads_file"]))
CONTIGS_FILES  = dict(zip(data["sample"], data["contigs_file"]))
KEGG_FILES     = dict(zip(data["sample"], data["kegg_file"]))

MODE_TO_FILENAME = {
            "pathway":   "pathway_abundance",
            "organisms": "organisms_abundance",
            "taxonomy":  "taxonomy",
}
# Inverse pour récupérer le mode depuis le wildcard filename
FILENAME_TO_MODE = {v: k for k, v in MODE_TO_FILENAME.items()}

READS_TREATMENT = config["output_path"]["treatment"] + "reads/"
CONTIGS_TREATMENT = config["output_path"]["treatment"] + "contigs/"
KEGG_TREATMENT = config["output_path"]["treatment"] + "kegg/"

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
# ==========================================================================
# Helpers
# ==========================================================================
def filter_active_sources(sources_list):
    """Filtre la liste des sources (reads, contigs, kegg) selon les run_xxx actifs"""
    filtered = []
    for s in sources_list:
        if s == "reads" and config.get("run_reads", False):
            filtered.append(s)
        elif s == "contigs" and config.get("run_contigs", False):
            filtered.append(s)
        elif s == "kegg" and config.get("run_kegg", False):
            filtered.append(s)
    return filtered

def phyloseq(pattern, sources_key):
    if not config["run_phyloseq"]:
        return []
    return expand(pattern, source=config["datatypes"][sources_key])

def deseq2(pattern, sources_key):
    if not config["run_deseq2"]:
        return []
    return expand(pattern, source=config["datatypes"][sources_key])

# ==========================================================================
# Cibles finales
# ==========================================================================
def get_targets():
    targets = []

    # --- Taxonomy ---
    if config["run_taxonomy"]:
        targets.extend([
            config["input_path"]["taxonomy_ncbi"]["local_path"], 
            config["input_path"]["pathway_bakta"]["local_path"],
        ])

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

    # --- Plots ---
    if config["run_plots"]:
        # Heatmap
        targets += phyloseq(
            config["output_path"]["plots"] + "heatmap/Heatmap_{source}.pdf",
            "heatmap"
        )
            
        # Stackedbarplots standard
        for mode, sources in config["datatypes"]["stackedbarplot"]["standard"].items():
            active_sources = filter_active_sources(sources)
            targets += expand(
                config["output_path"]["plots"] + f"stackedbarplot/{mode}/Stackedbarplot_{MODE_TO_FILENAME[mode]}_{{source}}.pdf",
                source=active_sources
            )

        # Stackedbarplots DESeq2
        targets += deseq2(
            config["output_path"]["plots"] + "stackedbarplot/deseq2/Stackedbarplot_deseq2_{source}.pdf",
            "deseq2"
        )

        # Volcano
        targets += deseq2(
            config["output_path"]["plots"] + "volcano/{source}/Volcano_deseq2_{source}.pdf",
            "volcano"
        )

    # --- Phyloseq ---
    if config["run_phyloseq"]:
        datatypes = ["reads", "contigs", "kegg"]
        active_phyloseq_dt = filter_active_sources(datatypes)
        
        for d in active_phyloseq_dt:
            targets.append(
            config["output_path"]["rds"] + f"{d}/phyloseq_{d}.rds"
            )

    # --- DESeq2 ---
    if config["run_deseq2"]:
        datatypes = ["contigs", "kegg"]
        active_deseq_dt = filter_active_sources(datatypes)
        contrasts = config.get("deseq2_contrasts", ["ref", "date", "combo"])
        for d in active_deseq_dt:
            for c in contrasts:
                targets.append(
                    config["output_path"]["rds"] + f"{d}/deseq2_{c}_{d}.rds"
                )
                targets.append(
                    config["output_path"]["parquet"] + f"{d}/deseq2_{c}_{d}.parquet"
                )

    if config["run_physico"]:
        targets.append(
            config["output_path"]["plots"] + "physico/Physico_plots.pdf"
        )
    # --- QC ---
    if config["run_qc"]:
        datatypes = ["reads", "contigs"]  # KEGG non intégré pour le moment
        active_qc_dt = filter_active_sources(datatypes)
        
        if active_qc_dt:  # au moins une source active
            targets.extend(expand([
                config["output_path"]["qc"] + "qc/Report_QC_final_{source}.pdf",
                config["output_path"]["parquet"] + "report_qc_final_{source}.parquet"
                ],
                source=active_qc_dt)
            )

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
        taxonomy = config["input_path"]["taxonomy_ncbi"]["local_path"]
    params:
        url      = config["input_path"]["taxonomy_ncbi"]["zip_url"]
    conda:   "../envs/py_env.yaml"
    script: "../scripts/utils/utils_convert_NCBInames_to_TaxaTable.py"

rule download_input_pathway:
    output:
        tsv = config["input_path"]["pathway_bakta"]["local_path"]
    params:
        url = config["input_path"]["pathway_bakta"]["url"]
    conda:   "../envs/py_env.yaml"
    script: "../scripts/utils/utils_pathway_levels_extraction.py"

# ==========================================================================
# UTILS — QC
# ==========================================================================
rule run_qc:
    input:
        data = lambda w: TREATMENT_SOURCES[w.source]
    output:
        pdf     = config["output_path"]["qc"] + "qc/Report_QC_final_{source}.pdf",
        parquet = config["output_path"]["parquet"] + "report_qc_final_{source}.parquet"
    params:
        active_modules = lambda w: filter_active_sources(["reads", "contigs"])
    conda:  "../envs/r_env.yaml"
    script: "../scripts/utils/utils_qc_wrapper.py"

# ==========================================================================
# READS — 3 étapes
# ==========================================================================

rule reads_countig:
    input:
        raw_data    = lambda w: READS_FILES[w.sample]
    output:
        counted     = READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/reads/01_Reads_counting_raw.py"

rule reads_filter:
    input:
        data     = READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    output:
        filtered = READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv"
    conda:   "../envs/py_env.yaml"
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
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/reads/03_Reads_add_taxaname.py"

# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================

rule contigs_counting:
    input:
        raw_data            = lambda w: CONTIGS_FILES[w.sample]
    output:
        counted             = CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    params:
        length_threshold    = config["contigs"]["length_threshold"]
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/contigs/01_Contigs_counting.py"

rule contigs_filter:
    input:
        data                = CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
        #current_sample = "results/contigs/1.Counted/counted_{sample}_contigs.tsv"
    output:
        filtered            = CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv"
    params:
        abundance_threshold = config["contigs"]["abundance_threshold"]
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/contigs/02_Contigs_filter.py"

rule contigs_rpkm:
    input:
        data = CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv"
    output:
        rpkm = CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/contigs/03_Contigs_RPKM.py"

rule contigs_rpkm_filter:
    input:
        data            = CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv"
    output:
        rpkm_filtered   = CONTIGS_TREATMENT + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    params:
        rpkm_threshold  = config["contigs"]["rpkm_threshold"]
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/contigs/04_Contigs_RPKM_filter.py"

rule contigs_union:
    input:
        abundance       = CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv",
        rpkm_filtered   = CONTIGS_TREATMENT + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    output:
        union           = CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/contigs/05_Contigs_union_filtered.py"

rule contigs_add_taxaname:
    input:
        data        = CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv",
        taxonomy    = config["input_path"]["taxonomy_megahit"]
    output:
        taxaname    = CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/contigs/06_Contigs_add_taxaname.py"

# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================

rule kegg_extraction:
    input:
        raw_data    = lambda w: KEGG_FILES[w.sample]
    output:
        extracted   = KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/kegg/01_Kegg_extraction.py"

rule kegg_intersec_count:
    input:
        data        = KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv",
        counted     = CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    output:
        intersec    = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/kegg/02_Kegg_count_intersection.py"

rule kegg_prepared_deseq2:
    input:
        data    = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        deseq2  = KEGG_TREATMENT + "3.Deseq2/deseq2_{sample}_kegg.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/kegg/03_Kegg_prepared_for_DESeq2.py"

rule kegg_standardization:
    input:
        data            = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv",
        current_sample  = KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        stand           = KEGG_TREATMENT + "021.Standardized/standardized_{sample}_kegg.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/kegg/021_Kegg_standardization.py"

rule kegg_standardization_agregation:
    input:
        data    = KEGG_TREATMENT + "021.Standardized/standardized_{sample}_kegg.tsv"
    output:
        agreg   = KEGG_TREATMENT + "022.Aggregated/stand_aggreg_{sample}_kegg.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/kegg/022_Kegg_standardization_aggregation.py"

rule kegg_merge_input_pathway_levels:
    input:
        data        = KEGG_TREATMENT + "022.Aggregated/stand_aggreg_{sample}_kegg.tsv",
        input_pathway     = config["input_path"]["pathway_bakta"]["local_path"]
    output:
        taxaname    = KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv"
    conda:   "../envs/py_env.yaml"
    script:  "../scripts/kegg/023_Kegg_merge_pathway_levels.py"

# ==========================================================================
# PLOTS R
# ==========================================================================

rule plot_stackedbarplot_deseq2:
    input:
        deseq_files = lambda w: expand(
            config["output_path"]["rds"] + "{source}/deseq2_{contrast}_{source}.rds",
            source=w.source,
            contrast=config.get("deseq2_contrasts", ["ref", "date", "combo"])
        ),
        phyloseq    = config["output_path"]["rds"] + "{source}/phyloseq_{source}.rds",
        metadata    = config["input_path"]["metadata"]
    output:
        pdf     = config["output_path"]["plots"] + "stackedbarplot/deseq2/Stackedbarplot_deseq2_{source}.pdf",
        parquet = config["output_path"]["parquet"] + "stackedbarplot_deseq2_{source}.parquet"
    params:
        padj      = config["plots"]["volcano"]["pvalue_threshold"],
        lfc       = config["plots"]["volcano"]["lfc_treshold"],
        contrasts = config.get("deseq2_contrasts", ["ref", "date", "combo"])
    conda:  "../envs/r_env.yaml"
    script: "../scripts/plots/Stackedbarplot_from_DESeq2.R"

rule plot_stackedbarplot:
    input:
        data     = lambda w: TREATMENT_SOURCES[w.source],
        metadata = config["input_path"]["metadata"]
    output:
        pdf     = config["output_path"]["plots"] + "stackedbarplot/{mode}/Stackedbarplot_{filename}_{source}.pdf",
        parquet = config["output_path"]["parquet"] + "stackedbarplot/{mode}/stackedbarplot_{filename}_{source}.parquet"
    wildcard_constraints:
        # Empêche les wildcards {mode} et {filename} de capturer le mot "deseq2"
        mode = "(?!deseq2)[a-zA-Z0-9_]+",
        filename = "(?!deseq2)[a-zA-Z0-9_]+"
    params:
        mode        = lambda w: FILENAME_TO_MODE.get(w.filename, w.filename),
        top_n       = config["plots"]["stackedbarplot"]["top_n"],
        target_rank = config["plots"]["stackedbarplot"]["taxon_rank"],
        value_col   = lambda w: PLOT_PARAMS[w.source]["value_col"]
    conda:  "../envs/r_env.yaml"
    script: "../scripts/plots/Stackedbarplot_unified.R"

rule plot_heatmap:
    input:
        rds             = config["output_path"]["rds"] + "{source}/phyloseq_{source}.rds",
        data            = lambda w: TREATMENT_SOURCES[w.source],
        metadata        = config["input_path"]["metadata"]
    output:
        pdf             = config["output_path"]["plots"] + "heatmap/Heatmap_{source}.pdf",
        parquet         = config["output_path"]["parquet"] + "heatmap_{source}.parquet",
    params:
        shared          = config["plots"]["shared"],
        top_n           = config["plots"]["heatmap"]["top_n"],
        clust_method    = config["plots"]["heatmap"]["clust_method"]
    conda:  "../envs/r_env.yaml"
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
    conda:  "../envs/r_env.yaml"
    script: "../scripts/plots/PCA.R"

rule plot_volcano_DESeq2:
    input:
        deseq_files = lambda w: expand(
            config["output_path"]["rds"] + "{source}/deseq2_{contrast}_{source}.rds",
            source=[w.source],
            contrast=config.get("deseq2_contrasts", ["ref", "date", "combo"])
        )
    output:
        pdf         = config["output_path"]["plots"] + "volcano/{source}/Volcano_deseq2_{source}.pdf",
        parquet     = config["output_path"]["parquet"] + "{source}/volcano_from_deseq2_{source}.parquet"
    params:
        padj        = config["plots"]["volcano"]["pvalue_threshold"],
        lfc         = config["plots"]["volcano"]["lfc_treshold"]
    conda:  "../envs/r_env.yaml"
    script: "../scripts/plots/Volcano_from_DESeq2.R"

rule plot_physico:
    input:
        physico  = config["input_path"]["physico_params"]
    output:
        pdf     = config["output_path"]["plots"] + "physico/Physico_plots.pdf",
        parquet = config["output_path"]["parquet"] + "physico/physico_data.parquet"
    conda:  "../envs/r_env.yaml"
    script: "../scripts/plots/Curves_physico_parameters.R"

# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================

rule run_phyloseq:
    input:
        data        = lambda w: TREATMENT_SOURCES[w.d],
        metadata    = config["input_path"]["metadata"]
    output:
        rds         = config["output_path"]["rds"] + "{d}/phyloseq_{d}.rds",
    conda:  "../envs/r_env.yaml"
    script: "../scripts/analysis/Phyloseq.R"

rule run_deseq2:
    input:
        data        = lambda w: TREATMENT_SOURCES["kegg_deseq2" if w.d == "kegg" else w.d],
        metadata    = config["input_path"]["metadata"]
    output:
        rds         = config["output_path"]["rds"] + "{d}/deseq2_{c}_{d}.rds",
        parquet     = config["output_path"]["parquet"] + "{d}/deseq2_{c}_{d}.parquet"

    conda:  "../envs/r_env.yaml"
    script: "../scripts/analysis/DESeq2.R"