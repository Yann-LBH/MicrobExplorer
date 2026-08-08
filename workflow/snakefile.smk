import os
import csv
#import pandas as pd
from pathlib import Path
from glob import glob


# Dev sous Windows patch os.path.join
def pjoin(*args):

    return os.path.normpath(os.path.join(*args))


configfile: "config/config.yaml"

# ==========================================================================
#   DICTIONARIES
# ==========================================================================

# Parse sample IDs without importing pandas at parsing time
with open(config["samples_file"], "r") as f:
    reader = csv.reader(f, delimiter="\t")
    next(reader)  # Skip header
    SAMPLES = [row[0] for row in reader if row]
#DATA = pd.read_table(config["samples_file"], index_col=0) #get the samples ID from the first column of the config file
#SAMPLES = DATA.index.tolist()
READS_FILES = {sample: f"{config['input_path']['data_raw']['reads']}reads_{sample}.kaijuNR" for sample in SAMPLES}
CONTIGS_FILES = {sample: f"{config['input_path']['data_raw']['contigs']}count-contigs-coassembly-{sample}.tsv" for sample in SAMPLES}
KEGG_FILES = f"{config['input_path']['data_raw']['kegg']}coassembly_bakta.gff3"

READS_TREATMENT = config["output_path"]["treatment"] + "reads/"
CONTIGS_TREATMENT = config["output_path"]["treatment"] + "contigs/"
KEGG_TREATMENT = config["output_path"]["treatment"] + "kegg/"

TAXONOMY = {
    "reads":config["input_path"]["taxonomy_ncbi"]["local_path"],
    "contigs":config["input_path"]["taxonomy_megahit"]["converted"],
    "kegg":config["input_path"]["pathway_bakta"]["local_path"]
}

TREATMENT_SOURCES = {
    "reads": expand(
        READS_TREATMENT + "3.CPM/cpm_{sample}_reads.tsv", sample=SAMPLES
    ),
    "reads_deseq": expand(
        READS_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_reads.tsv", sample=SAMPLES
    ),
    "reads_phyloseq": expand(
        READS_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_reads.tsv", sample=SAMPLES
    ),
    "contigs": expand(
        CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "contigs_deseq": expand(
        CONTIGS_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "contigs_phyloseq": expand(
        CONTIGS_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "kegg": expand(
        KEGG_TREATMENT + "3.RPKM/rpkm_{sample}_kegg.tsv", sample=SAMPLES
    ),
    "kegg_deseq": expand(
        KEGG_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_kegg.tsv", sample=SAMPLES,
    ),
    "kegg_phyloseq": expand(
        KEGG_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_kegg.tsv", sample=SAMPLES,
    )
}

PARQUET_FILES = {
    "stackedbarplot": glob(pjoin(config["output_path"]["parquet"], "stackedbarplot", "*.parquet")),
    "heatmap": glob(pjoin(config["output_path"]["parquet"], "heatmap", "*.parquet")),
    "pca": glob(pjoin(config["output_path"]["parquet"], "pca", "*.parquet")),
    "volcano": glob(pjoin(config["output_path"]["parquet"], "volcano", "*.parquet")),
    "physico": glob(pjoin(config["output_path"]["parquet"], "physico", "*.parquet")),
    "qc": glob(pjoin(config["output_path"]["parquet"], "qc", "*.parquet"))
}

QC_STEPS = {
"reads": {
        # Raw data path
        "brut": (pjoin(config["input_path"]["data_raw"]["reads"]), "reads_{sample}.kaijuNR"),
        # Treatment output paths
        "counted": (pjoin(config["output_path"]["treatment"], "reads", "1.Counted"), "counted_{sample}_reads.tsv"),
        "filtered": (pjoin(config["output_path"]["treatment"], "reads", "2.Filtered"), "filtered_{sample}_reads.tsv"),
        "cpm": (pjoin(config["output_path"]["treatment"], "reads", "3.CPM"), "cpm_{sample}_reads.tsv")
    },
    "contigs": {
        # Raw data path
        "brut": (pjoin(config["input_path"]["data_raw"]["contigs"]), "count-contigs-coassembly-{sample}.tsv"),
        # Treatment output paths
        "counted": (pjoin(config["output_path"]["treatment"], "contigs", "1.Counted"), "counted_{sample}_contigs.tsv"),
        "filtered": (pjoin(config["output_path"]["treatment"], "contigs", "2.Filtered"), "filtered_{sample}_contigs.tsv"),
        "rpkm": (pjoin(config["output_path"]["treatment"], "contigs", "3.RPKM"), "rpkm_{sample}_contigs.tsv"),
        "rpkm_filtered": (pjoin(config["output_path"]["treatment"], "contigs", "4.RPKM_Filtered"), "rpkm_filtered_{sample}_contigs.tsv"),
        "union": (pjoin(config["output_path"]["treatment"], "contigs", "5.Union"), "union_{sample}_contigs.tsv"),
        "matrix_deseq": (pjoin(config["output_path"]["treatment"], "contigs", "Matrix", "Deseq"), "matrix_deseq_{sample}_contigs.tsv"),
        "matrix_phyloseq": (pjoin(config["output_path"]["treatment"], "contigs", "Matrix", "Phyloseq"), "matrix_phyloseq_{sample}_contigs.tsv")
    },
    "kegg": {
        # Raw annotations or input files
        "brut": (pjoin(config["input_path"]["data_raw"]["kegg"]), "coassembly_bakta.gff3"),
        # Treatment and filtering steps
        "extracted": (pjoin(config["output_path"]["treatment"], "kegg", "1.Extracted"), "extracted_{sample}_kegg.tsv"),
        "intersected": (pjoin(config["output_path"]["treatment"], "kegg", "2.Intersected"), "intersected_{sample}_kegg.tsv"),
        "rpkm": (pjoin(config["output_path"]["treatment"], "kegg", "3.RPKM"), "rpkm_{sample}_kegg.tsv"),
        "standardized": (pjoin(config["output_path"]["treatment"], "kegg", "4.Standardized"), "standardized_{sample}_kegg.tsv"),
        "aggregated": (pjoin(config["output_path"]["treatment"], "kegg", "5.Aggregated"), "stand_aggreg_{sample}_kegg.tsv"),
        "matrix_deseq": (pjoin(config["output_path"]["treatment"], "kegg", "Matrix", "Deseq"), "matrix_deseq_{sample}_kegg.tsv"),
        "matrix_phyloseq": (pjoin(config["output_path"]["treatment"], "kegg", "Matrix", "Phyloseq"), "matrix_phyloseq_{sample}_kegg.tsv")
    },
}

PLOT_PARAMS = {
    "reads": {"stand_col": "cpm", "label": "cpm"},
    "contigs": {"stand_col": "rpkm", "label": "rpkm"},
    "kegg": {"stand_col": "rpkm", "label": "rpkm"},
}


# ==========================================================================
# Helpers
# ==========================================================================
def filter_active_sources(sources_list):
    """Filters the list of data sources (reads, contigs, kegg) based on active run_xxx flags."""

    return [source for source in sources_list if config.get(f"run_{source}", False)]


def phyloseq(pattern, sources_key):
    if not config["run_phyloseq"]:
        return []
    active_sources = filter_active_sources(config["datatypes"][sources_key])
    return expand(pattern, source=active_sources)


def deseq2(pattern, sources_key):
    if not config["run_deseq2"]:
        return []
    active_sources = filter_active_sources(config["datatypes"][sources_key])
    return expand(pattern, source=active_sources)


def get_graphs_input(wildcards):
    return TREATMENT_SOURCES[wildcards.source]


def pca(pattern, sources_key):
    """Generates targets for PCA if the run_pca flag is active in the config."""
    if not config.get("run_pca", True):  # Default to True if the flag isn't set yet
        return []
    active_sources = filter_active_sources(config["datatypes"][sources_key])
    return expand(pattern, source=active_sources)


def get_qc_inputs(source):
    input_files = []
    
    for step, (folder, fname) in QC_STEPS[source].items():
        # If the filename is static (like coassembly_bakta.gff3), add it only once
        if "{sample}" not in fname:
            input_files.append(pjoin(folder, fname))
        else:
            # Otherwise, expand it for all samples
            for s in SAMPLES:
                input_files.append(pjoin(folder, fname.format(sample=s)))
                
    return input_files

# ==========================================================================
# Cibles finales
# ==========================================================================
def get_targets():
    targets = []

    # --- Taxonomy ---
    if config.get("run_taxonomy", False):
        targets.extend([
                config["input_path"]["taxonomy_ncbi"]["local_path"],
                config["input_path"]["taxonomy_megahit"]["converted"],
                config["input_path"]["pathway_bakta"]["local_path"],
        ])

    # --- Reads ---
    if config.get("run_reads", False):
        targets.extend(TREATMENT_SOURCES["reads"])

    # --- Contigs ---
    if config.get("run_contigs", False):
        targets.extend(TREATMENT_SOURCES["contigs"])
        targets.extend(TREATMENT_SOURCES["contigs_deseq"])
        targets.extend(TREATMENT_SOURCES["contigs_phyloseq"])

    if config.get("run_kegg", False):
        targets.extend(TREATMENT_SOURCES["kegg"])
        targets.extend(TREATMENT_SOURCES["kegg_deseq"])
        targets.extend(TREATMENT_SOURCES["kegg_phyloseq"])

    # --- Plots ---
    if config.get("run_plots", False):
        # Stackedbarplots standard
        for mode, sources in config["datatypes"]["stackedbarplot"]["standard"].items():
            active_sources = filter_active_sources(sources)
            targets.extend(
                expand(
                    [
                        pjoin(
                            config["output_path"]["plots"],
                            "stackedbarplot",
                            "{mode}",
                            "Stackedbarplot_{mode}_{source}.pdf",
                        ),
                        pjoin(
                            config["output_path"]["parquet"],
                            "stackedbarplot",
                            "stackedbarplot_{mode}_{source}.parquet",
                        ),
                    ],
                        mode=mode,
                        source=active_sources,
                )
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
                "stackedbarplot",
                "stackedbarplot_deseq2_{source}.parquet",
            ),
            "deseq2",
        )

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
            pjoin(config["output_path"]["parquet"], "heatmap","heatmap_{source}.parquet"),
            "heatmap",
        )

        # PCA
        targets += pca(
            pjoin(config["output_path"]["plots"], "pca", "PCA_{source}.pdf"),
            "pca"
        )
        targets += pca(
            pjoin(config["output_path"]["parquet"], "pca", "pca_{source}.parquet"),
            "pca"
        )
        targets += pca(
            pjoin(config["output_path"]["plots"], "pca", "pca_contributions_{source}.xlsx"),
            "pca"
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
                "volcano",
                "volcano_from_deseq2_{source}.parquet",
            ),
            "volcano",
        )

    # --- DESeq2 ---
    if config.get("run_deseq2", False):
        targets += deseq2(
            pjoin(config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"),
            "deseq2",
        )
        targets += deseq2(
            pjoin(
                config["output_path"]["parquet"], "deseq2", "deseq2_{source}.parquet"
            ),
            "deseq2",
        )

    # --- Phyloseq ---
    if config.get("run_phyloseq", False):
        targets += phyloseq(
            pjoin(config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds"),
            "phyloseq",
        )

    if config.get("run_physico", False):
        targets.extend(
            expand(
                [
                    pjoin(config["output_path"]["plots"], "physico", "Physico_plot.pdf"),
                    pjoin(config["output_path"]["parquet"], "physico", "physico_plot.parquet"),
                ]
            )
        )
    # --- QC ---
    if config.get("run_qc", False):
        datatypes = list(QC_STEPS.keys())
        active_qc_dt = filter_active_sources(datatypes)
        active_qc_dt = [
            source for source in active_qc_dt
            if len(get_qc_inputs(source)) > 0
        ]

        if active_qc_dt:
            targets.extend(
                expand(
                    [
                        pjoin(
                            config["output_path"]["qc"],
                            "Report_QC_final_{source}.pdf",
                        ),
                        pjoin(
                            config["output_path"]["parquet"],
                            "qc",
                            "report_qc_final_{source}.parquet",
                        ),
                    ],
                    source=active_qc_dt,
                )
            )

#    if config.get("run_shiny", False):       
#        targets.extend(
#            expand(
#                [
#                    pjoin(config["output_path"]["parquet"], "shiny_stackedbarplot.parquet"),
#                    pjoin(config["output_path"]["parquet"], "shiny_heatmap.parquet"),
#                    pjoin(config["output_path"]["parquet"], "shiny_pca.parquet"),
#                    pjoin(config["output_path"]["parquet"], "shiny_volcano.parquet"),
#                    pjoin(config["output_path"]["parquet"], "shiny_physico.parquet"),
#                    pjoin(config["output_path"]["parquet"], "shiny_qc.parquet")
#                ]
#            )
#        )

rule all:
    input:
        #config["output_path"]["audit"],
        "benchmarks/summary_benchmarks.csv"
        get_targets(),
# ==========================================================================
# UTILS — Taxonomy and input_pathways
# ==========================================================================
rule download_taxonomy:
    output:
        zip=temp("data/taxonomy/new_taxdump.zip"),
        taxonomy=TAXONOMY["reads"]
    conda:
        "envs/py_env.yaml"
    params:
        url=config["input_path"]["taxonomy_ncbi"]["zip_url"],
        dmp_name=config["input_path"]["taxonomy_ncbi"]["dmp_name"],
    script:
        os.path.abspath(
            "workflow/scripts/utils/utils_convert_NCBInames_to_TaxaTable.py"
        )

rule convert_megahit_taxonomy:
    input:
        megahit=config["input_path"]["taxonomy_megahit"]["original"]
    output:
        taxonomy=TAXONOMY["contigs"]
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath(
            "workflow/scripts/utils/utils_format_MegahitTable.py"
        )

rule download_pathway:
    output:
        pathway=TAXONOMY["kegg"]
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
            config["output_path"]["parquet"], "qc", "report_qc_data_{source}.parquet"
        )
    conda:
        "envs/py_env.yaml"
    params:
        steps_config=lambda w: QC_STEPS[w.source],
        active_modules=filter_active_sources(config["datatypes"]["qc"])
    script:
        os.path.abspath("workflow/scripts/utils/utils_qc_wrapper.py")


rule run_plot_qc:
    input:
        data=pjoin(config["output_path"]["parquet"], "qc", "report_qc_data_{source}.parquet")
    output:
        pdf=pjoin(config["output_path"]["qc"], "Report_QC_final_{source}.pdf"),
        parquet=pjoin(config["output_path"]["parquet"], "qc", "report_qc_final_{source}.parquet")
    conda:
        "envs/r_env.yaml"
    params:
        steps_config=lambda w: QC_STEPS[w.source],
        active_modules=filter_active_sources(config["datatypes"]["qc"])
    script:
        os.path.abspath("workflow/scripts/utils/utils_barplot_qc.R")


# ==========================================================================
# READS — 3 étapes
# ==========================================================================


rule reads_counting:
    input:
        raw_data=lambda w: READS_FILES[w.sample]
    output:
        counted=READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    benchmark:
        "benchmarks/reads/{sample}_counted.tsv"
    #threads: 4
    #resources:
        #mem_mb=8000,       # 8 Go de RAM
        #runtime=60         # 60 minutes max (utile pour SLURM/SGE)
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/reads/01_Reads_counting.py")


rule reads_filter:
    input:
        data=READS_TREATMENT + "1.Counted/counted_{sample}_reads.tsv"
    output:
        filtered=READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv"
    benchmark:
        "benchmarks/reads/{sample}_filtered.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        count_threshold=config["reads"]["count_threshold"]
    script:
        os.path.abspath("workflow/scripts/reads/02_Reads_filter.py")


rule reads_CPM:
    input:
        data=READS_TREATMENT + "2.Filtered/filtered_{sample}_reads.tsv"
    output:
        cpm=READS_TREATMENT + "3.CPM/cpm_{sample}_reads.tsv"
    benchmark:
        "benchmarks/reads/{sample}_cpm.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/reads/03_Reads_CPM.py")


# ==========================================================================
# CONTIGS — 6 étapes
# ==========================================================================


rule contigs_counting:
    input:
        raw_data=lambda w: CONTIGS_FILES[w.sample]
    output:
        counted=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    benchmark:
        "benchmarks/contigs/{sample}_counted.tsv"
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
    benchmark:
        "benchmarks/contigs/{sample}_global_abundance.tsv"
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
    benchmark:
        "benchmarks/contigs/{sample}_filtered.tsv"
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
    benchmark:
        "benchmarks/contigs/{sample}_rpkm.tsv"
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
    benchmark:
        "benchmarks/contigs/{sample}_rpkm_filtered.tsv"
    conda:
        "envs/py_env.yaml"
    params:
        rpkm_threshold=config["contigs"]["rpkm_threshold"]
    script:
        os.path.abspath("workflow/scripts/contigs/04_Contigs_RPKM_filter.py")


rule contigs_union:
    input:
        data=CONTIGS_TREATMENT + "3.RPKM/rpkm_{sample}_contigs.tsv",
        abundance=CONTIGS_TREATMENT + "2.Filtered/filtered_{sample}_contigs.tsv",
        rpkm_filtered=CONTIGS_TREATMENT
        + "4.RPKM_Filtered/rpkm_filtered_{sample}_contigs.tsv"
    output:
        union=CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv"
    benchmark:
        "benchmarks/contigs/{sample}_union.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/contigs/05_Contigs_union_filtered.py")


rule contigs_create_matrix:
    input:
        data=CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv"
    output:
        matrix_deseq=CONTIGS_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_contigs.tsv",
        matrix_phyloseq=CONTIGS_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_contigs.tsv"
    benchmark:
        "benchmarks/contigs/{sample}_create_matrix.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/contigs/06_Contigs_create_matrix.py")


# ==========================================================================
# KEGG — 6 étapes
# ==========================================================================


rule kegg_extraction:
    input:
        raw_data=KEGG_FILES
    output:
        extracted=KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv"
    benchmark:
        "benchmarks/kegg/{sample}_extraction.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/01_Kegg_extraction.py")


rule kegg_intersec_count:
    input:
        data=KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv",
        counted=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    output:
        intersec=KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv",
        matrix_deseq=KEGG_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_kegg.tsv"
    benchmark:
        "benchmarks/kegg/{sample}_intersec_count.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/02_Kegg_count_intersec_and_aggregate.py")


rule kegg_rpkm:
    input:
        data=KEGG_TREATMENT + "2.Intersected/intersected_{sample}_kegg.tsv"
    output:
        rpkm=KEGG_TREATMENT + "3.RPKM/rpkm_{sample}_kegg.tsv"
        matrix_phyloseq=KEGG_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_kegg.tsv"
    benchmark:
        "benchmarks/kegg/{sample}_rpkm.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/03_Kegg_RPKM_and_aggregate.py")


# ==========================================================================
# PLOTS R
# ==========================================================================


rule plot_stackedbarplot:
    input:
        data=get_graphs_input,
        metadata=config["input_path"]["metadata"],
        taxonomy=lambda w: TAXONOMY[w.source],
        title_resolver = config["title_resolver"]
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "stackedbarplot",
            "{mode}",
            "Stackedbarplot_{mode}_{source}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "stackedbarplot",
            "stackedbarplot_{mode}_{source}.parquet",
        )
    wildcard_constraints:
        # Prevents the {mode} and {source} wildcards from matching the word “deseq2” and “kegg_stand”
        mode = "(?!deseq2)[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        language=config["language"],
        shared=config["plots"]["shared"],
        mode=lambda w: w.mode,
        top_n=config["plots"]["stackedbarplot"]["top_n"],
        stand_col=lambda w: PLOT_PARAMS[w.source]["stand_col"],
        rank=lambda w: config["plots"]["shared"]["rank"][w.source],
    script:
        os.path.abspath("workflow/scripts/plots/Stackedbarplot_abundance.R")


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
            config["output_path"]["parquet"], 
            "stackedbarplot",   
            "stackedbarplot_deseq2_{source}.parquet"
        )
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title=lambda w: config["plots"]["stackedbarplot_deseq2"]["title"],
        subtitle=lambda w: config["plots"]["stackedbarplot_deseq2"]["subtitle"],
        contrast=config["deseq2"]["contrast"],
        padj=config["plots"]["stackedbarplot_deseq2"]["contrasts_values"]["padj_threshold"],
        lfc=config["plots"]["stackedbarplot_deseq2"]["contrasts_values"]["lfc_threshold"],
        top_n=config["plots"]["stackedbarplot_deseq2"]["top_n"],
        rank=lambda w: config["plots"]["shared"]["rank"][w.source],
    script:
        os.path.abspath("workflow/scripts/plots/Stackedbarplot_from_DESeq2.R")


rule plot_heatmap:
    input:
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds"
        ),
        title_resolver = config["title_resolver"]
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "heatmap",
            "{source}",
            "Heatmap_{source}.pdf",
        ),
        parquet=pjoin(config["output_path"]["parquet"], "heatmap", "heatmap_{source}.parquet")
    conda:
        "envs/r_env.yaml"
    params:
        language=config["language"],
        shared=config["plots"]["shared"],
        top_n=config["plots"]["heatmap"]["top_n"],
        clust_method=config["plots"]["heatmap"]["clust_method"],
        distance_method=config["plots"]["heatmap"]["distance_method"],
        rank=lambda w: config["plots"]["shared"]["rank"][w.source],
    script:
        os.path.abspath("workflow/scripts/plots/Heatmap.R")


rule plot_pca:
    input:
        data=get_graphs_input,
        metadata=config["input_path"]["metadata"],
        physico=config["input_path"]["physico_params"],
        taxonomy=lambda w: TAXONOMY[w.source],
        title_resolver = config["title_resolver"]
    output:
        pdf=pjoin(config["output_path"]["plots"], "pca", "PCA_{source}.pdf"),
        parquet=pjoin(config["output_path"]["parquet"], "pca", "pca_{source}.parquet"),
        xlsx=pjoin(config["output_path"]["plots"], "pca", "pca_contributions_{source}.xlsx")
    conda:
        "envs/r_env.yaml"
    params:
        language=config["language"],
        shared=config["plots"]["shared"],
        title=lambda w: config["plots"]["pca"]["title"],
        subtitle=lambda w: config["plots"]["pca"]["subtitle"],
        top_n=config["plots"]["pca"]["top_n"],
        dim_x=config["plots"]["pca"]["dim_x"],
        dim_y=config["plots"]["pca"]["dim_y"],
        physico_col=config["plots"]["pca"]["physico_col"],
        stand_col=lambda w: PLOT_PARAMS[w.source]["stand_col"],
        point_size=config["plots"]["pca"]["point_size"],
        rank=lambda w: config["plots"]["shared"]["rank"][w.source],
    script:
        os.path.abspath("workflow/scripts/plots/PCA.R")


rule plot_volcano_DESeq2:
    input:
        deseq_files=pjoin(
            config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"
        ),
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds"
        ),
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "volcano",
            "{source}",
            "Volcano_deseq2_{source}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "volcano",
            "volcano_from_deseq2_{source}.parquet",
        )
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title=lambda w: config["plots"]["volcano"]["title"],
        subtitle=lambda w: config["plots"]["volcano"]["subtitle"],
        contrast=config["deseq2"]["contrast"],
        padj=config["plots"]["stackedbarplot_deseq2"]["contrasts_values"]["padj_threshold"],
        lfc=config["plots"]["stackedbarplot_deseq2"]["contrasts_values"]["lfc_threshold"],
        top_n=config["plots"]["volcano"]["top_n"],
        rank=lambda w: config["plots"]["shared"]["rank"][w.source],
    script:
        os.path.abspath("workflow/scripts/plots/Volcano_from_DESeq2.R")


rule plot_physico:
    input:
        physico_params=config["input_path"]["physico_params"]
    output:
        pdf=pjoin(config["output_path"]["plots"], "physico", "Physico_plot.pdf"),
        parquet=pjoin(
            config["output_path"]["parquet"], "physico", "physico_plot.parquet"
        )
    conda:
        "envs/r_env.yaml"
    script:
        os.path.abspath("workflow/scripts/plots/Curves_physico_parameters.R")


# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================

rule run_deseq2:
    input:
        data=lambda w: TREATMENT_SOURCES[f"{w.source}_deseq"],
        metadata=config["input_path"]["metadata"],
        taxonomy=lambda w: TAXONOMY[w.source],
    output:
        rds=pjoin(config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"),
        parquet=pjoin(config["output_path"]["parquet"], "deseq2", "deseq2_{source}.parquet")
    conda:
        "envs/r_env.yaml"
    params:
        contrast=config["deseq2"]["contrast"],
        ref=config["deseq2"]["ref"],
        sizefactor=config["deseq2"]["sizefactor"],
        test=config["deseq2"]["test"],
        fittype=config["deseq2"]["fittype"],
    script:
        os.path.abspath("workflow/scripts/analysis/DESeq2.R")



rule run_phyloseq:
    input:
        data=lambda w: TREATMENT_SOURCES[f"{w.source}_phyloseq"],
        metadata=config["input_path"]["metadata"],
        taxonomy=lambda w: TAXONOMY[w.source],
    output:
        rds=pjoin(config["output_path"]["rds"], "{source}", "phyloseq_{source}.rds")
    params:
        stand_col=lambda w: PLOT_PARAMS[w.source]["stand_col"],
    conda:
        "envs/r_env.yaml"
    script:
        os.path.abspath("workflow/scripts/analysis/Phyloseq.R")


# ==========================================================================
# SHINY MASTER PARQUET by Analyse
# ==========================================================================

#rule shiny_master_parquet:
#    input:
#        stackedbarplot = PARQUET_FILES["stackedbarplot"],
#        heatmap        = PARQUET_FILES["heatmap"],
#        pca            = PARQUET_FILES["pca"],
#        volcano        = PARQUET_FILES["volcano"],
#        physico        = PARQUET_FILES["physico"],
#        qc             = PARQUET_FILES["qc"]
#    output:
#        stackedbarplot = pjoin(config["output_path"]["parquet"], "shiny_stackedbarplot.parquet"),
#        heatmap        = pjoin(config["output_path"]["parquet"], "shiny_heatmap.parquet"),
#        pca            = pjoin(config["output_path"]["parquet"], "shiny_pca.parquet"),
#        volcano        = pjoin(config["output_path"]["parquet"], "shiny_volcano.parquet"),
#        physico        = pjoin(config["output_path"]["parquet"], "shiny_physico.parquet"),
#        qc             = pjoin(config["output_path"]["parquet"], "shiny_qc.parquet")
#    conda:
#        "envs/py_env.yaml"
#    script:
#        os.path.abspath("workflow/scripts/utils/utils_shiny_master_parquet.py")

# ==========================================================================
# Benchmarks stats
# ==========================================================================


rule aggregate_benchmarks:
    input:
        reads=expand("benchmarks/reads/{sample}_counting.tsv", sample=SAMPLES),
        contigs=expand("benchmarks/contigs/{sample}_counting.tsv", sample=SAMPLES),
        kegg=expand("benchmarks/kegg/{sample}_extraction.tsv", sample=SAMPLES)
    output:
        csv="benchmarks/summary_benchmarks.csv"
    run:
        import pandas as pd
        from pathlib import Path

        records = []
        for file_path in input:
            path = Path(file_path)
            df = pd.read_csv(path, sep="\t")
            
            # Extract metadata from file structure
            df["rule"] = path.parent.name
            df["sample"] = path.stem.replace("_counting", "")
            records.append(df)

        # Combine all benchmark data into a single dataframe
        summary_df = pd.concat(records, ignore_index=True)
        
        # Reorder key columns first
        primary_cols = ["rule", "sample", "s", "h:m:s", "max_rss", "cpu_time"]
        other_cols = [c for c in summary_df.columns if c not in primary_cols]
        summary_df[primary_cols + other_cols].to_csv(output.csv, index=False)
# ==========================================================================
# AUDIT
# ==========================================================================

ACTIVE_SOURCES = filter_active_sources(["reads", "contigs", "kegg"])
include: "audit/audit.smk"