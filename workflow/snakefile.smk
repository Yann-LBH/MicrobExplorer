import os
import pandas as pd
from pathlib import Path


# Dev sous Windows patch os.path.join
def pjoin(*args):

    return os.path.normpath(os.path.join(*args))


configfile: "config/config.yaml"


# ==========================================================================
#   DICTIONARIES
# ==========================================================================

# Échantillons
DATA = pd.read_table(config["samples_file"], index_col=0) #get the samples ID from the first column of the config file
SAMPLES = DATA.index.tolist()
READS_FILES = {sample: f"{config['input_path']['data_raw']['reads']}reads_{sample}.kaijuNR" for sample in SAMPLES}
CONTIGS_FILES = {sample: f"{config['input_path']['data_raw']['contigs']}count-contigs-coassembly-{sample}.tsv" for sample in SAMPLES}
KEGG_FILES = {sample: f"{config['input_path']['data_raw']['kegg']}kegg_{sample}.gff3" for sample in SAMPLES}

MODE_TO_FILENAME = {
    "pathway": "pathway_abundance",
    "organisms": "organisms_abundance",
    "taxonomy": "taxonomy",
}
# Inverse pour récupérer le mode depuis le wildcard filename
FILENAME_TO_MODE = {v: k for k, v in MODE_TO_FILENAME.items()}

READS_TREATMENT = config["output_path"]["treatment"] + "reads/"
CONTIGS_TREATMENT = config["output_path"]["treatment"] + "contigs/"
KEGG_TREATMENT = config["output_path"]["treatment"] + "kegg/"

TREATMENT_SOURCES = {
    "reads": expand(
        READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv", sample=SAMPLES
    ),
    "contigs": expand(
        CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "kegg": expand(
        KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv", sample=SAMPLES
    ),
    "kegg_deseq2": expand(
        KEGG_TREATMENT + "3.Deseq2/deseq2_{sample}_kegg.tsv", sample=SAMPLES
    ),
}

QC_STEPS = {
    "reads": {
        "brut": ("Data", "reads_{sample}.kaijuNR"),
        "counted": ("1.Counted", "counted_{sample}_reads.tsv"),
        "filtered": ("2.Filtered", "filtered_{sample}_reads.tsv"),
        "final": ("3.Annotated", "annotated_{sample}_reads.tsv"),
    },
    "contigs": {
        "brut": ("Data", "count-contigs-coassembly-{sample}.tsv"), #need adaptation
        "counted": ("1.Counted", "counted_{sample}_contigs.tsv"),
        "filtered": ("2.Filtered", "filtered_{sample}_contigs.tsv"),
        "rpkm": ("3.Rpkm", "rpkm_{sample}_contigs.tsv"),
        "rpkm_filtered": ("4.Rpkm_filtered", "rpkm_filtered_{sample}_contigs.tsv"),
        "union": ("5.Union", "union_{sample}_contigs.tsv"),
        "final": ("6.Annotated", "annotated_{sample}_contigs.tsv"),
    },
}

PLOT_PARAMS = {
    "reads": {"value_col": "Reads", "label": "Read counts"},
    "contigs": {"value_col": "RPKM", "label": "RPKM"},
    "kegg": {"value_col": "RPKM", "label": "RPKM"},
}


# ==========================================================================
# Helpers
# ==========================================================================
def filter_active_sources(sources_list):
    """Filters the list of data sources (reads, contigs, kegg) based on active run_xxx flags."""
    # Une seule ligne de compréhension de liste remplace tout le bloc if/elif
    return [source for source in sources_list if config.get(f"run_{source}", False)]


def phyloseq(pattern, sources_key):
    if not config["run_phyloseq"]:
        return []
    return expand(pattern, source=config["datatypes"][sources_key])


def deseq2(pattern, sources_key):
    if not config["run_deseq2"]:
        return []
    active_sources = filter_active_sources(config["datatypes"][sources_key])
    return expand(pattern, source=active_sources)


def get_qc_inputs(source):
    return [
        pjoin(folder, fname.format(sample=s))
        for step, (folder, fname) in QC_STEPS[source].items()
        for s in SAMPLES
    ]


# ==========================================================================
# Cibles finales
# ==========================================================================
def get_targets():
    targets = []

    # --- Taxonomy ---
    if config["run_taxonomy"]:
        targets.extend(
            [
                config["input_path"]["taxonomy_ncbi"]["local_path"],
                config["input_path"]["pathway_bakta"]["local_path"],
            ]
        )

    # --- Reads ---
    if config["run_reads"]:
        targets.extend(
            expand(
                READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv",
                sample=SAMPLES,
            )
        )

    # --- Contigs ---
    if config["run_contigs"]:
        targets.extend(
            expand(
                CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv",
                sample=SAMPLES,
            )
        )

    if config["run_kegg"]:
        targets.extend(
            expand(
                KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv",
                sample=SAMPLES,
            )
        )

    # --- Plots ---
    if config["run_plots"]:
        # Heatmap
        targets += phyloseq(
            pjoin(
                config["output_path"]["plots"],
                "heatmap",
                "{source}",
                "Heatmap_{source}.pdf",
            ),
            "heatmap",
        )
        targets += phyloseq(
            pjoin(config["output_path"]["parquet"], "heatmap_{source}.parquet"),
            "heatmap",
        )

        # Stackedbarplots standard
        for mode, sources in config["datatypes"]["stackedbarplot"]["standard"].items():
            active_sources = filter_active_sources(sources)
            targets += expand(
                pjoin(
                    config["output_path"]["plots"],
                    "stackedbarplot",
                    mode,
                    "Stackedbarplot_" + MODE_TO_FILENAME[mode] + "_{source}.pdf",
                ),
                source=active_sources,
            )

        # Stackedbarplots DESeq2
        targets += deseq2(
            pjoin(
                config["output_path"]["plots"],
                "stackedbarplot",
                "deseq2",
                "Stackedbarplot_deseq2_{source}.pdf",
            ),
            "deseq2",
        )
        targets += deseq2(
            pjoin(
                config["output_path"]["parquet"],
                "stackedbarplot_deseq2_{source}.parquet",
            ),
            "deseq2",
        )

        # Volcano
        targets += deseq2(
            pjoin(
                config["output_path"]["plots"],
                "volcano",
                "{source}",
                "Volcano_deseq2_{source}.pdf",
            ),
            "volcano",
        )
        targets += deseq2(
            pjoin(
                config["output_path"]["parquet"],
                "{source}",
                "volcano_from_deseq2_{source}.parquet",
            ),
            "volcano",
        )

    # --- Phyloseq ---
    if config["run_phyloseq"]:
        targets += phyloseq(
            pjoin(config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds"),
            "phyloseq",
        )

    # --- DESeq2 ---
    if config["run_deseq2"]:
        targets += deseq2(
            pjoin(config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"),
            "deseq2",
        )
        targets += deseq2(
            pjoin(
                config["output_path"]["parquet"], "{source}", "deseq2_{source}.parquet"
            ),
            "deseq2",
        )

    if config["run_physico"]:
        targets.append(
            pjoin(config["output_path"]["plots"], "physico", "Physico_plots.pdf")
        )
    # --- QC ---
    if config["run_qc"]:
        datatypes = list(
            QC_STEPS.keys()
        )  # ["reads", "contigs"] — tiré du dict, pas en dur

        active_qc_dt = [
            source
            for source in datatypes
            if source in QC_STEPS  # étape définie
            and config.get("sources", {}).get(source, False)  # activée dans config
            and len(get_qc_inputs(source)) > 0  # fichiers réellement attendus
        ]

        if active_qc_dt:
            targets.extend(
                expand(
                    [
                        pjoin(
                            config["output_path"]["qc"],
                            "qc",
                            "Report_QC_final_{source}.pdf",
                        ),
                        pjoin(
                            config["output_path"]["parquet"],
                            "report_qc_final_{source}.parquet",
                        ),
                    ],
                    source=active_qc_dt,
                )
            )

    return targets


rule all:
    input:
        "logs/final_pipeline_report.txt",
        get_targets(),


# ==========================================================================
# UTILS — Taxonomy and input_pathways
# ==========================================================================
rule download_taxonomy:
    output:
        zip=temp("data/taxonomy/new_taxdump.zip"),
        taxonomy=config["input_path"]["taxonomy_ncbi"]["local_path"]
    conda:
        "envs/py_env.yaml"
    params:
        url=config["input_path"]["taxonomy_ncbi"]["zip_url"],
                dmp_name=config["input_path"]["taxonomy_ncbi"]["dmp_name"],
    script:
        os.path.abspath(
            "workflow/scripts/utils/utils_convert_NCBInames_to_TaxaTable.py"
        )


rule download_input_pathway:
    output:
        tsv=config["input_path"]["pathway_bakta"]["local_path"]
    conda:
        "envs/py_env.yaml"
    params:
        url=config["input_path"]["pathway_bakta"]["url"]
    script:
        os.path.abspath("workflow/scripts/utils/utils_pathway_levels_extraction.py")


# ==========================================================================
# UTILS — QC
# ==========================================================================
rule run_qc:
    input:
        data=lambda w: get_qc_inputs(w.source)
    output:
        parquet=pjoin(
            config["output_path"]["parquet"], "report_qc_data_{source}.parquet"
        )
    conda:
        "envs/py_env.yaml"
    params:
        steps_config=lambda w: QC_STEPS[w.source],
        active_modules=lambda w: filter_active_sources(config["datatypes"]["qc"])
    script:
        os.path.abspath("workflow/scripts/utils/utils_qc_wrapper.py")


rule run_plot_qc:
    input:
        data=pjoin(config["output_path"]["parquet"], "report_qc_data_{source}.parquet")
    output:
        pdf=pjoin(config["output_path"]["qc"], "qc", "Report_QC_final_{source}.pdf"),
        parquet=pjoin(
            config["output_path"]["parquet"], "report_qc_final_{source}.parquet"
        )
    conda:
        "envs/r_env.yaml"
    params:
        active_modules=lambda w: filter_active_sources(config["datatypes"]["qc"])
    script:
        os.path.abspath("workflow/scripts/utils/utils_barplot_qc.R")


# ==========================================================================
# READS — 3 étapes
# ==========================================================================


rule reads_countig:
    input:
        raw_data=lambda w: READS_FILES[w.sample]
    output:
        counted=READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/reads/01_Reads_counting.py")


rule reads_filter:
    input:
        data=READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    output:
        filtered=READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        count_threshold=config["reads"]["count_threshold"]
    script:
        os.path.abspath("workflow/scripts/reads/02_Reads_filter.py")


rule reads_add_taxaname:
    input:
        data=READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv",
        taxonomy=config["input_path"]["taxonomy_ncbi"]["local_path"]
    output:
        taxaname=READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        top_n=config["reads"]["top_n"],
        noise=config["reads"]["noise"]
    script:
        os.path.abspath("workflow/scripts/reads/03_Reads_add_taxaname.py")


# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================


rule contigs_counting:
    input:
        raw_data=lambda w: CONTIGS_FILES[w.sample]
    output:
        counted=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        length_threshold=config["contigs"]["length_threshold"]
    script:
        os.path.abspath("workflow/scripts/contigs/01_Contigs_counting.py")


# Optimisation of rule filter use for global abundance calculation
rule contigs_global_abundance:
    input:
        all_data=expand(
            CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv",
            sample=SAMPLES,
        )
    output:
        global_abundance=CONTIGS_TREATMENT + "2.Global_Abundance/global_abundance_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/contigs/02_a_Contigs_global_abundance.py")


rule contigs_filter:
    input:
        data=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv",
        global_abundance=CONTIGS_TREATMENT + "2.Global_Abundance/global_abundance_contigs.tsv"
    output:
        filtered=CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        abundance_threshold=config["contigs"]["abundance_threshold"]
    script:
        os.path.abspath("workflow/scripts/contigs/02_b_Contigs_filter.py")


rule contigs_rpkm:
    input:
        data=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    output:
        rpkm=CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/contigs/03_Contigs_RPKM.py")


rule contigs_rpkm_filter:
    input:
        data=CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv"
    output:
        rpkm_filtered=CONTIGS_TREATMENT
        + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        rpkm_threshold=config["contigs"]["rpkm_threshold"],
    script:
        os.path.abspath("workflow/scripts/contigs/04_Contigs_RPKM_filter.py")


rule contigs_union:
    input:
        raw_data=lambda w: CONTIGS_FILES[w.sample],
        abundance=CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv",
        rpkm_filtered=CONTIGS_TREATMENT
        + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    output:
        union=CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/contigs/05_Contigs_union_filtered.py")


rule contigs_add_taxaname:
    input:
        data=CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv",
        taxonomy=config["input_path"]["taxonomy_megahit"]
    output:
        taxaname=CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/contigs/06_Contigs_add_taxaname.py")


# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================


rule kegg_extraction:
    input:
        raw_data=lambda w: KEGG_FILES[w.sample]
    output:
        extracted=KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/01_Kegg_extraction.py")


rule kegg_intersec_count:
    input:
        data=KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv",
        counted=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    output:
        intersec=KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/02_Kegg_count_intersection.py")


rule kegg_prepared_deseq2:
    input:
        data=KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        deseq2=KEGG_TREATMENT + "3.Deseq2/deseq2_{sample}_kegg.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/03_Kegg_prepared_for_DESeq2.py")


rule kegg_standardization:
    input:
        data=KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        stand=KEGG_TREATMENT + "021.Standardized/standardized_{sample}_kegg.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/021_Kegg_standardization.py")


rule kegg_standardization_agregation:
    input:
        data=KEGG_TREATMENT + "021.Standardized/standardized_{sample}_kegg.tsv"
    output:
        agreg=KEGG_TREATMENT + "022.Aggregated/stand_aggreg_{sample}_kegg.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/022_Kegg_standardization_aggregation.py")


rule kegg_merge_input_pathway_levels:
    input:
        data=KEGG_TREATMENT + "022.Aggregated/stand_aggreg_{sample}_kegg.tsv",
        input_pathway=config["input_path"]["pathway_bakta"]["local_path"]
    output:
        taxaname=KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/023_Kegg_merge_pathway_levels.py")


# ==========================================================================
# PLOTS R
# ==========================================================================


rule plot_stackedbarplot_deseq2:
    input:
        deseq_files=pjoin(
            config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"
        ),
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds"
        ),
        metadata=config["input_path"]["metadata"]
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "stackedbarplot",
            "deseq2",
            "Stackedbarplot_deseq2_{source}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"], "stackedbarplot_deseq2_{source}.parquet"
        )
    conda:
        "envs/r_env.yaml"
    params:
        title=config["plots"]["stackedbarplot"]["title"],
        subtitle=config["plots"]["stackedbarplot"]["subtitle"],
        padj=config["plots"]["volcano"]["pvalue_threshold"],
        lfc=config["plots"]["volcano"]["lfc_treshold"],
        contrasts=config["deseq2_contrasts"],
        top_n=config["plots"]["stackedbarplot"]["top_n"]
    script:
        os.path.abspath("workflow/scripts/plots/Stackedbarplot_from_DESeq2.R")


rule plot_stackedbarplot:
    input:
        data=lambda w: TREATMENT_SOURCES[w.source],
        metadata=config["input_path"]["metadata"]
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "stackedbarplot",
            "{mode}",
            "Stackedbarplot_{filename}_{source}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "stackedbarplot",
            "{mode}",
            "stackedbarplot_{filename}_{source}.parquet",
        )
    wildcard_constraints:
        # Empêche les wildcards {mode} et {filename} de capturer le mot "deseq2"
        mode="(?!deseq2)[a-zA-Z0-9_]+",
        filename="(?!deseq2)[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        mode=lambda w: FILENAME_TO_MODE.get(w.filename, w.filename),
        top_n=config["plots"]["stackedbarplot"]["top_n"],
        target_rank=config["plots"]["stackedbarplot"]["taxon_rank"],
        value_col=lambda w: PLOT_PARAMS[w.source]["value_col"]
    script:
        os.path.abspath("workflow/scripts/plots/Stackedbarplot_abundance.R")


rule plot_heatmap:
    input:
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds"
        ),
        data=lambda w: TREATMENT_SOURCES[w.source],
        metadata=config["input_path"]["metadata"]
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "heatmap",
            "{source}",
            "Heatmap_{source}.pdf",
        ),
        parquet=pjoin(config["output_path"]["parquet"], "heatmap_{source}.parquet")
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        top_n=config["plots"]["heatmap"]["top_n"],
        clust_method=config["plots"]["heatmap"]["clust_method"]
    script:
        os.path.abspath("workflow/scripts/plots/Heatmap.R")


rule plot_pca:
    input:
        data=lambda w: TREATMENT_SOURCES[w.source],
        metadata=config["input_path"]["metadata"],
        physico=config["input_path"]["physico_params"]
    output:
        pdf=pjoin(config["output_path"]["plots"], "pca", "PCA_{source}.pdf"),
        parquet=pjoin(config["output_path"]["parquet"], "pca_{source}.parquet"),
        csv=pjoin("results", "infos", "pca", "pca_contributions_{source}.csv")
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        top_n=config["plots"]["pca"]["top_n"],
        point_size=config["plots"]["pca"]["point_size"],
        dim_x=config["plots"]["pca"]["dim_x"],
        dim_y=config["plots"]["pca"]["dim_y"],
        physico_cols=config["plots"]["pca"]["physico_cols"]
    script:
        os.path.abspath("workflow/scripts/plots/PCA.R")


rule plot_volcano_DESeq2:
    input:
        deseq_files=pjoin(
            config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"
        )
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "volcano",
            "{source}",
            "Volcano_deseq2_{source}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "{source}",
            "volcano_from_deseq2_{source}.parquet",
        )
    conda:
        "envs/r_env.yaml"
    params:
        padj=config["plots"]["volcano"]["pvalue_threshold"],
        lfc=config["plots"]["volcano"]["lfc_treshold"],
        contrast=config["deseq2_contrasts"]
    script:
        os.path.abspath("workflow/scripts/plots/Volcano_from_DESeq2.R")


rule plot_physico:
    input:
        physico_params=config["input_path"]["physico_params"]
    output:
        pdf=pjoin(config["output_path"]["plots"], "physico", "Physico_plots.pdf"),
        parquet=pjoin(
            config["output_path"]["parquet"], "physico", "physico_data.parquet"
        )
    conda:
        "envs/r_env.yaml"
    script:
        os.path.abspath("workflow/scripts/plots/Curves_physico_parameters.R")


# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================


rule run_phyloseq:
    input:
        data=lambda w: TREATMENT_SOURCES[w.source],
        metadata=config["input_path"]["metadata"]
    output:
        rds=pjoin(config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds")
    conda:
        "envs/r_env.yaml"
    script:
        os.path.abspath("workflow/scripts/analysis/Phyloseq.R")


rule run_deseq2:
    input:
        data=lambda w: TREATMENT_SOURCES[
            "kegg_deseq2" if w.source == "kegg" else w.source
        ],
        metadata=config["input_path"]["metadata"]
    output:
        rds=pjoin(config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"),
        parquet=pjoin(
            config["output_path"]["parquet"], "{source}", "deseq2_{source}.parquet"
        )
    conda:
        "envs/r_env.yaml"
    params:
        contrasts=config["deseq2_contrasts"]
    script:
        os.path.abspath("workflow/scripts/analysis/DESeq2.R")


# ==========================================================================
# FINAL REPORT - PIPELINE AUDIT
# ==========================================================================


rule generate_pipeline_report:
    input:
        # Enlèvement de la virgule de fin pour éviter le tuple fantôme
        get_targets()
    output:
        report="logs/final_pipeline_report.txt"
    run:
        import os

        # English comments in script blocks
        failed_files = {}

        def check_files(category, file_list):
            """Helper to audit a list of paths and log missing ones."""
            missing = [
                f
                for f in file_list
                if not os.path.exists(f) or os.path.getsize(f) == 0
            ]
            if missing:
                failed_files[category] = missing

        # CORRECTION : Utilisation de snakemake.config au lieu de config
        # --- 1. Taxonomy ---
        if snakemake.config["run_taxonomy"]:
            check_files(
                "Taxonomy Databases",
                [
                    snakemake.config["input_path"]["taxonomy_ncbi"]["local_path"],
                    snakemake.config["input_path"]["pathway_bakta"]["local_path"],
                ],
            )
        # --- 2. Reads ---
        if snakemake.config["run_reads"]:
            check_files(
                "Reads Annotations (TSV)",
                expand(
                    READS_TREATMENT + "3.Annotated/annotated_{sample}_reads.tsv",
                    sample=SAMPLES,
                ),
            )
        # --- 3. Contigs ---
        if snakemake.config["run_contigs"]:
            check_files(
                "Contigs Annotations (TSV)",
                expand(
                    CONTIGS_TREATMENT
                    + "6.Annotated/annotated_{sample}_contigs.tsv",
                    sample=SAMPLES,
                ),
            )
        # --- 4. KEGG ---
        if snakemake.config["run_kegg"]:
            check_files(
                "KEGG Annotations (TSV)",
                expand(
                    KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv",
                    sample=SAMPLES,
                ),
            )
        # --- 5. Phyloseq Objects ---
        if snakemake.config["run_phyloseq"]:
            check_files(
                "Phyloseq (RDS)",
                phyloseq(
                    pjoin(
                        snakemake.config["output_path"]["rds"],
                        "{source}",
                        "phyloseq_{source}.rds",
                    ),
                    "phyloseq",
                ),
            )
        # --- 6. DESeq2 Objects & Parquets ---
        if snakemake.config["run_deseq2"]:
            check_files(
                "DESeq2 (RDS/Parquet)",
                deseq2(
                    pjoin(
                        snakemake.config["output_path"]["rds"],
                        "{source}",
                        "deseq2_{source}.rds",
                    ),
                    "deseq2",
                )
                + deseq2(
                    pjoin(
                        snakemake.config["output_path"]["parquet"],
                        "{source}",
                        "deseq2_{source}.parquet",
                    ),
                    "deseq2",
                ),
            )
        # --- 7. Plots & Associated Parquets ---
        if snakemake.config["run_plots"]:
            plot_targets = []
            # Heatmaps
            plot_targets += phyloseq(
                pjoin(
                    snakemake.config["output_path"]["plots"],
                    "heatmap",
                    "{source}",
                    "Heatmap_{source}.pdf",
                ),
                "heatmap",
            )
            plot_targets += phyloseq(
                pjoin(snakemake.config["output_path"]["parquet"], "heatmap_{source}.parquet"),
                "heatmap",
            )
            # Stackedbarplots standard
            for mode, sources in snakemake.config["datatypes"]["stackedbarplot"][
                "standard"
            ].items():
                active_sources = filter_active_sources(sources)
                plot_targets += expand(
                    pjoin(
                        snakemake.config["output_path"]["plots"],
                        "stackedbarplot",
                        mode,
                        "Stackedbarplot_"
                        + MODE_TO_FILENAME[mode]
                        + "_{source}.pdf",
                    ),
                    source=active_sources,
                )
            # Stackedbarplots DESeq2
            plot_targets += deseq2(
                pjoin(
                    snakemake.config["output_path"]["plots"],
                    "stackedbarplot",
                    "deseq2",
                    "Stackedbarplot_deseq2_{source}.pdf",
                ),
                "deseq2",
            )
            plot_targets += deseq2(
                pjoin(
                    snakemake.config["output_path"]["parquet"],
                    "stackedbarplot_deseq2_{source}.parquet",
                ),
                "deseq2",
            )
            # Volcano
            plot_targets += deseq2(
                pjoin(
                    snakemake.config["output_path"]["plots"],
                    "volcano",
                    "{source}",
                    "Volcano_deseq2_{source}.pdf",
                ),
                "volcano",
            )
            plot_targets += deseq2(
                pjoin(
                    snakemake.config["output_path"]["parquet"],
                    "{source}",
                    "volcano_from_deseq2_{source}.parquet",
                ),
                "volcano",
            )

            check_files("Plots & Figure Parquets", plot_targets)
        # --- 8. Physico ---
        if snakemake.config["run_physico"]:
            check_files(
                "Physico Plots",
                [
                    pjoin(
                        snakemake.config["output_path"]["plots"],
                        "physico",
                        "Physico_plots.pdf",
                    )
                ],
            )
        # --- 9. Quality Control (QC) ---
        if snakemake.config["run_qc"]:
            datatypes = list(QC_STEPS.keys())
            active_qc_dt = [
                s
                for s in datatypes
                if s in QC_STEPS
                and snakemake.config.get("sources", {}).get(s, False)
                and len(get_qc_inputs(s)) > 0
            ]
            if active_qc_dt:
                check_files(
                    "QC Reports & Parquets",
                    expand(
                        [
                            pjoin(
                                snakemake.config["output_path"]["qc"],
                                "qc",
                                "Report_QC_final_{source}.pdf",
                            ),
                            pjoin(
                                snakemake.config["output_path"]["parquet"],
                                "report_qc_final_{source}.parquet",
                            ),
                        ],
                        source=active_qc_dt,
                    ),
                )
        
        # CORRECTION : open(snakemake.output.report) au lieu de open(output.report)
        # --- Write the final markdown-like report ---
        with open(snakemake.output.report, "w") as f:
            f.write("==================================================\n")
            f.write("         PIPELINE RUN AUDIT REPORT                \n")
            f.write("==================================================\n\n")

            if not failed_files:
                f.write(
                    "🎉 SUCCESS: All requested steps completed with zero errors!\n"
                )
                f.write(
                    "All expected samples, matrix parquets, RDS objects, and PDFs are valid.\n"
                )
            else:
                f.write(
                    "⚠️ WARNING: The pipeline finished with missing or failed targets.\n"
                )
                f.write(
                    "Because some samples failed, final global steps (like plots) may have been skipped.\n\n"
                )
                f.write("--- DETAILED MISSING TARGETS PER CATEGORY ---\n")
                for category, files in failed_files.items():
                    f.write(f"\n❌ {category} ({len(files)} missing file(s)):\n")
                    for file in files:
                        f.write(f"  - {file}\n")
                        
        # Terminal feedback
        if failed_files:
            print(
                f"\n⚠️ [Snakemake Audit] Some steps failed. Check the summary here: {snakemake.output.report}\n"
            )
        else:
            print(
                f"\n🎉 [Snakemake Audit] Perfect run! Summary generated: {snakemake.output.report}\n"
            )