import os
import csv
import re
import json
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
    "reads_counts": expand(
        READS_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_reads.tsv", sample=SAMPLES
    ),
    "reads_normalized": expand(
        READS_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_reads.tsv", sample=SAMPLES
    ),
    "contigs": expand(
        CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "contigs_counts": expand(
        CONTIGS_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "contigs_normalized": expand(
        CONTIGS_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_contigs.tsv", sample=SAMPLES
    ),
    "kegg": expand(
        KEGG_TREATMENT + "3.RPKM/rpkm_{sample}_kegg.tsv", sample=SAMPLES
    ),
    "kegg_counts": expand(
        KEGG_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_kegg.tsv", sample=SAMPLES,
    ),
    "kegg_normalized": expand(
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
        "cpm": (pjoin(config["output_path"]["treatment"], "reads", "3.CPM"), "cpm_{sample}_reads.tsv"),
        "matrix_deseq": (pjoin(config["output_path"]["treatment"], "reads", "Matrix", "Deseq"), "matrix_deseq_{sample}_reads.tsv"),
        "matrix_phyloseq": (pjoin(config["output_path"]["treatment"], "reads", "Matrix", "Phyloseq"), "matrix_phyloseq_{sample}_reads.tsv")
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

METRIC_PARAMS = {
    "counts": {
        "reads":   {"stand_col": "count",       "label": "Counts"},
        "contigs": {"stand_col": "read_mapped", "label": "Mapped Reads"},
        "kegg":    {"stand_col": "read_mapped", "label": "Mapped Reads"},
    },
    "normalized": {
        "reads":   {"stand_col": "cpm",  "label": "CPM"},
        "contigs": {"stand_col": "rpkm", "label": "RPKM"},
        "kegg":    {"stand_col": "rpkm", "label": "RPKM"},
    }
}

KEGG_RANKS = ["level_3", "gene_description"]

RULE_DATATYPE_MAP = {
    "phyloseq": "phyloseq",
    "deseq2": "deseq2",
    "stackedbarplot": "stackedbarplot",
    "heatmap": "heatmap",
    "volcano": "volcano",
    "pca": "pca",
}

BASE_SOURCES = ["reads", "contigs", "kegg"]

# ==========================================================================
# Helpers
# ==========================================================================
# Multilingual traduction Title and Subtitle | I18N
with open(config.get("title_resolver"), "r", encoding="utf-8") as f:
    translations_data = json.load(f)["translation"]
    I18N = {item["key"]: item for item in translations_data}

def get_text(key):
    """Retrieves the template string in the active language."""
    lang = config.get("language", "en")
    return I18N.get(key, {}).get(lang, "")

def is_enabled(flag_name):
    """Utility function to check if a module flag is active"""
    return config.get(f"run_{flag_name}", False)

def filter_active_sources(sources_list):
    """Filters the list of data sources (reads, contigs, kegg) based on active run_xxx flags."""
    return [source for source in sources_list if is_enabled(source)]

def get_ranks_for_source(source, plot_type=None):
    ranks = config["plots"]["shared"]["rank"].get(source, [])
    if isinstance(ranks, str):
        ranks = [ranks]

    if source == "kegg":
        if plot_type == "heatmap_deseq2":
            # La heatmap prend level_3 ET gene_description
            return ranks
        else:
            # Tous les autres plots (PCA, Volcano, etc.) excluent level_3
            return [r for r in ranks if r != "level_3"]

    # Pour reads/contigs, retourne tous les rangs taxonomiques (species, genus...)
    return ranks

def is_clr_mode(rule_name="permanova"):
    """Vérifie si le mode CLR est activé pour une règle donnée."""
    cfg = config.get(rule_name, {})
    use_clr = cfg.get("use_clr", False)
    dist_method = str(cfg.get("distance_method", "")).lower()
    return use_clr or dist_method == "clr"

def get_permanova_inputs(wildcards):
    """Sélectionne dynamiquement les tables (counts vs normalized)."""
    metric = "counts" if is_clr_mode("permanova") else "normalized"
    suffix = "_counts" if metric == "counts" else "_normalized"
    return TREATMENT_SOURCES[f"{wildcards.source}{suffix}"]

def get_permanova_params(wildcards):
    """Extrait la colonne d'abondance et le label associés."""
    metric = "counts" if is_clr_mode("permanova") else "normalized"
    return METRIC_PARAMS[metric][wildcards.source]

# --- Module Helpers ---

def taxonomy():
    if not is_enabled("taxonomy"):
        return []
    return [
        config["input_path"]["taxonomy_ncbi"]["local_path"],
        config["input_path"]["taxonomy_megahit"]["converted"],
        config["input_path"]["pathway_bakta"]["local_path"],
    ]

def reads():
    if not is_enabled("reads"):
        return []
    return TREATMENT_SOURCES["reads"] + TREATMENT_SOURCES["reads_counts"] + TREATMENT_SOURCES["reads_normalized"]

def contigs():
    if not is_enabled("contigs"):
        return []
    return TREATMENT_SOURCES["contigs"] + TREATMENT_SOURCES["contigs_counts"] + TREATMENT_SOURCES["contigs_normalized"]

def kegg():
    if not is_enabled("kegg"):
        return []
    return TREATMENT_SOURCES["kegg"] + TREATMENT_SOURCES["kegg_counts"] + TREATMENT_SOURCES["kegg_normalized"]

# --- Plots & Analysis Helpers ---

def physico():
    if not is_enabled("physico"):
        return []
    return [
        pjoin(config["output_path"]["plots"], "physico", "Physico_plot.pdf"),
        pjoin(config["output_path"]["parquet"], "physico", "physico_plot.parquet"),
    ]

def stackedbarplot_abundance():
    if not (is_enabled("stackedbarplot_abundance") and is_enabled("phyloseq")):
        return []

    targets = []
    abundance_dict = config.get("datatypes", {}).get("stackedbarplot", {}).get("abundance", {})

    # Parcourt chaque mode (pathway_relative, relative_by_sample, etc.) et filtre ses sources
    for mode, sources in abundance_dict.items():
        active_sources = filter_active_sources(sources)  # Correction : filtre la liste 'sources' du mode courant
        for source in active_sources:
            source_ranks = get_ranks_for_source(source, plot_type="stackedbarplot_abundance")
            targets.extend(
                expand(
                    [
                        pjoin(
                            config["output_path"]["plots"],
                            "stackedbarplot_abundance",
                            "{mode}",
                            "Stackedbarplot_{mode}_{source}_{rank}.pdf",
                        ),
                        pjoin(
                            config["output_path"]["parquet"],
                            "stackedbarplot_abundance",
                            "stackedbarplot_{mode}_{source}_{rank}.parquet",
                        ),
                    ],
                    mode=mode,
                    source=source,
                    rank=source_ranks,
                )
            )

    return targets


def stackedbarplot_deseq2():
    if not (is_enabled("stackedbarplot_deseq2") and is_enabled("phyloseq") and is_enabled("deseq2")):
        return []

    targets = []
    for source in filter_active_sources(config["datatypes"]["stackedbarplot"]["deseq2"]):
        # Taxonomic ranks are automatically retrieved for reads/contigs
        # and only gene_description for KEGG
        source_ranks = get_ranks_for_source(source, plot_type="stackedbarplot_deseq2")

        targets.extend(
            expand(
                [
                    pjoin(config["output_path"]["plots"], "stackedbarplot", "deseq2", "Stackedbarplot_deseq2_{source}_{rank}.pdf"),
                    pjoin(config["output_path"]["parquet"], "stackedbarplot", "stackedbarplot_deseq2_{source}_{rank}.parquet"),
                ],
                source=source,
                rank=source_ranks,
            )
        )

    return targets


def heatmap_abundance():
    # Correction : on vérifie 'run_heatmap' (conformément au config.yaml)
    if not (is_enabled("heatmap_abundance") and is_enabled("phyloseq")):
        return []

    targets = []
    sources = config.get("datatypes", {}).get("heatmap", {}).get("abundance", [])

    for source in filter_active_sources(sources):
        source_ranks = get_ranks_for_source(source, plot_type="heatmap_abundance")
        targets.extend(
            expand(
                [
                    pjoin(config["output_path"]["plots"], "heatmap", "{source}", "Heatmap_{source}_{rank}.pdf"),
                    pjoin(config["output_path"]["parquet"], "heatmap", "heatmap_{source}_{rank}.parquet"),
                ],
                source=source,
                rank=source_ranks,
            )
        )

    return targets


def heatmap_deseq2():
    if not (is_enabled("heatmap_deseq2") and is_enabled("phyloseq") and is_enabled("deseq2")):
        return []

    targets = []
    for source in filter_active_sources(config["datatypes"]["heatmap"]["deseq2"]):
        # Taxonomic ranks are automatically retrieved for reads/contigs
        # and only gene_description for KEGG
        source_ranks = get_ranks_for_source(source, plot_type="heatmap_deseq2")

        targets.extend(
            expand(
                [
                    pjoin(config["output_path"]["plots"], "heatmap", "{source}", "Heatmap_deseq2_{source}_{rank}.pdf"),
                    pjoin(config["output_path"]["parquet"], "heatmap", "heatmap_deseq2_{source}_{rank}.parquet"),
                ],
                source=source,
                rank=source_ranks,
            )
        )

    return targets


def pca():
    if not (is_enabled("pca") and is_enabled("phyloseq")):
        return []

    targets = []
    for source in filter_active_sources(config["datatypes"]["pca"]):
        # Taxonomic ranks are automatically retrieved for reads/contigs
        # and only gene_description for KEGG
        source_ranks = get_ranks_for_source(source, plot_type="pca")

        targets.extend(
            expand(
                [
                    pjoin(config["output_path"]["plots"], "pca", "{source}", "PCA_{source}_{rank}.pdf"),
                    pjoin(config["output_path"]["parquet"], "pca", "pca_{source}_{rank}.parquet"),
                    pjoin(config["output_path"]["plots"], "pca", "pca_contributions_{source}_{rank}.xlsx"),
                ],
                source=source,
                rank=source_ranks,
            )
        )

    return targets

def volcano():
    if not (is_enabled("volcano") and is_enabled("deseq2")):
        return []
        
    targets = []
    for source in filter_active_sources(config["datatypes"]["volcano"]):
        # Taxonomic ranks are automatically retrieved for reads/contigs
        # and only gene_description for KEGG
        source_ranks = get_ranks_for_source(source, plot_type="volcano")

        targets.extend(
            expand(
                [
                    pjoin(config["output_path"]["plots"], "volcano", "{source}", "Volcano_{source}_{rank}.pdf"),
                    pjoin(config["output_path"]["parquet"], "volcano", "volcano_{source}_{rank}.parquet"),
                ],
                source=source,
                rank=source_ranks,
            )
        )

    return targets

def permanova_analysis():
    if not is_enabled("permanova"):
        return []
        
    sources = config["datatypes"]["permanova"]
    active_sources = filter_active_sources(sources)
    if not active_sources:
        return []

    return expand(
        [
            pjoin(config["output_path"]["permanova"], "{source}", "Permanova_{source}.pdf"),
            pjoin(config["output_path"]["parquet"], "permanova", "permanova_{source}.parquet"),
        ],
        source=active_sources,
    )

def deseq2_analysis():
    if not is_enabled("deseq2"):
        return []
        
    sources = config["datatypes"]["deseq2"]
    active_sources = filter_active_sources(sources)
    if not active_sources:
        return []

    return expand(
        [
            pjoin(config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"),
            pjoin(config["output_path"]["parquet"], "deseq2", "deseq2_{source}.parquet"),
        ],
        source=active_sources,
    )

def phyloseq_analysis():
    if not is_enabled("phyloseq"):
        return []
        
    targets = []
    sources = config["datatypes"]["phyloseq"]
    active_sources = filter_active_sources(sources)

    for source in active_sources:
        ranks = config["plots"]["shared"]["rank"].get(source, [])
        if isinstance(ranks, str):
            ranks = [ranks]
        targets.extend(
            expand(
                pjoin(config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"),
                source=source,
                rank=ranks,
            )
        )

    return targets

def qc():
    if not is_enabled("qc"):
        return []
    datatypes = list(QC_STEPS.keys())
    active_qc_dt = [source for source in filter_active_sources(datatypes) if len(get_qc_inputs(source)) > 0]
    if not active_qc_dt:
        return []
    return expand(
        [
            pjoin(config["output_path"]["qc"], "Report_QC_final_{source}.pdf"),
            pjoin(config["output_path"]["parquet"], "qc", "report_qc_final_{source}.parquet"),
        ],
        source=active_qc_dt,
    )

def get_all_benchmarks(wildcards):
    if not is_enabled("benchmarks"):
        return []
    benchmarks = set()

    stacked_cfg = config.get("datatypes", {}).get("stackedbarplot", {})
    stacked_modes = dict(stacked_cfg.get("abundance", {}))
    if "deseq2" in stacked_cfg:
        stacked_modes["deseq2"] = stacked_cfg["deseq2"]

    for rule in workflow.rules:
        if not rule.benchmark:
            continue

        pattern = re.sub(r"\{(\w+)(?:,[^}]*)?\}", r"{\1}", str(rule.benchmark))
        names = set(re.findall(r"\{(\w+)\}", pattern))

        if not names:
            benchmarks.add(pattern)
            continue

        kwargs = {}

        # 1. Gestion de 'mode'
        if "mode" in names:
            kwargs["mode"] = list(stacked_modes.keys())

        # 2. Gestion de 'source'
        if "source" in names:
            datatype_key = RULE_DATATYPE_MAP.get(rule.name, "phyloseq")
            sources = config.get("datatypes", {}).get(datatype_key, BASE_SOURCES)
            active_sources = filter_active_sources(sources)
            if not active_sources:
                continue
            kwargs["source"] = active_sources

        # 3. Gestion de 'rank'
        if "rank" in names:
            sources_to_check = kwargs.get("source", BASE_SOURCES)
            all_ranks = set()
            for src in sources_to_check:
                r = config.get("plots", {}).get("shared", {}).get("rank", {}).get(src, ["all"])
                if isinstance(r, str):
                    all_ranks.add(r)
                else:
                    all_ranks.update(r)
            kwargs["rank"] = list(all_ranks)

        # 4. Gestion de 'sample'
        if "sample" in names:
            kwargs["sample"] = SAMPLES

        # S'assurer que tous les wildcards de la règle ont bien une valeur attribuée
        if names.issubset(kwargs.keys()):
            benchmarks.update(expand(pattern, **kwargs))

    return sorted(benchmarks)

# ==========================================================================
# Cibles finales
# ==========================================================================

def get_targets():
    targets = []

    # Collect targets sequentially from modular helper functions
    targets += taxonomy()
    targets += reads()
    targets += contigs()
    targets += kegg()

    # Plots
    targets += stackedbarplot_abundance()
    targets += stackedbarplot_deseq2()
    targets += heatmap_abundance()
    targets += heatmap_deseq2()
    targets += pca()
    targets += volcano()
    targets += physico()

    # Statistical Analyses & Data objects
    targets += permanova_analysis()
    targets += deseq2_analysis()
    targets += phyloseq_analysis()

    # Quality Control
    #targets += qc()

    return targets

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

#    return targets

rule all:
    input:
        #config["output_path"]["audit"],
        #"benchmarks/summary_benchmarks.csv",
        get_targets()
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
# READS — 4 étapes
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
        cpm=READS_TREATMENT + "3.CPM/cpm_{sample}_reads.tsv",
    benchmark:
        "benchmarks/reads/{sample}_cpm.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/reads/03_Reads_CPM.py")

rule reads_create_matrix:
    input:
        data=READS_TREATMENT + "3.CPM/cpm_{sample}_reads.tsv",
    output:
        matrix_deseq=READS_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_reads.tsv",
        matrix_phyloseq=READS_TREATMENT + "Matrix/Phyloseq/matrix_phyloseq_{sample}_reads.tsv"
    benchmark:
        "benchmarks/reads/{sample}_create_matrix.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/reads/04_Reads_create_matrix.py")


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
        "benchmarks/contigs/global_abundance.tsv"
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


rule kegg_union_count:
    input:
        data=KEGG_TREATMENT + "1.Extracted/extracted_{sample}_kegg.tsv",
        counted=CONTIGS_TREATMENT + "1.Counted/counted_{sample}_contigs.tsv"
    output:
        union=KEGG_TREATMENT + "2.Union/union_{sample}_kegg.tsv",
        matrix_deseq=KEGG_TREATMENT + "Matrix/Deseq/matrix_deseq_{sample}_kegg.tsv"
    benchmark:
        "benchmarks/kegg/{sample}_union_count.tsv"
    conda:
        "envs/py_env.yaml"
    script:
        os.path.abspath("workflow/scripts/kegg/02_Kegg_count_union_and_aggregate.py")


rule kegg_rpkm:
    input:
        data=KEGG_TREATMENT + "2.Union/union_{sample}_kegg.tsv"
    output:
        rpkm=KEGG_TREATMENT + "3.RPKM/rpkm_{sample}_kegg.tsv",
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


rule plot_stackedbarplot_abundance:
    input:
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"
        ),
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "stackedbarplot_abundance",
            "{mode}",
            "Stackedbarplot_{mode}_{source}_{rank}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "stackedbarplot_abundance",
            "stackedbarplot_{mode}_{source}_{rank}.parquet",
        ),
    benchmark:
        "benchmarks/{source}/{mode}_stackedbarplot_{rank}.tsv"
    wildcard_constraints:
        mode="|".join(config["datatypes"]["stackedbarplot"]["abundance"].keys()),
        source="reads|contigs|kegg",
        rank="[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title_template=lambda w: get_text("TITLE_STACKEDBARPLOT_ABUNDANCE"),
        subtitle_template=lambda w: get_text("SUBTITLE_STACKEDBARPLOT_ABUNDANCE"),
        text_scope_all=lambda w: get_text("all"),
        text_scope_each=lambda w: get_text("each"),
        mode=lambda w: w.mode,
        top_n=config["plots"]["stackedbarplot_abundance"]["top_n"],
        stand_col=lambda w: METRIC_PARAMS["normalized"][w.source]["stand_col"],
        rank=lambda w: w.rank,
    script:
        os.path.abspath("workflow/scripts/plots/Stackedbarplot_abundance.R")


rule plot_stackedbarplot_deseq2:
    input:
        deseq_files=pjoin(
            config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"
        ),
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"
        ),
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "stackedbarplot",
            "deseq2",
            "Stackedbarplot_deseq2_{source}_{rank}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "stackedbarplot",
            "stackedbarplot_deseq2_{source}_{rank}.parquet",
        ),
    benchmark:
        "benchmarks/{source}/stackedbarplot_deseq2_{rank}.tsv"
    wildcard_constraints:
        source="reads|contigs|kegg",
        rank="[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title_template=lambda w: get_text("TITLE_STACKEDBARPLOT_DESEQ2"),
        subtitle_template=lambda w: get_text("SUBTITLE_STACKEDBARPLOT_DESEQ2"),
        contrast=config["deseq2"]["contrast"],
        padj=config["plots"]["stackedbarplot_deseq2"]["contrasts_values"]["padj_threshold"],
        lfc=config["plots"]["stackedbarplot_deseq2"]["contrasts_values"]["lfc_threshold"],
        top_n=config["plots"]["stackedbarplot_deseq2"]["top_n"],
        rank=lambda w: w.rank,
    script:
        os.path.abspath("workflow/scripts/plots/Stackedbarplot_from_DESeq2.R")


rule plot_heatmap_abundance:
    input:
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"
        ),
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "heatmap",
            "{source}",
            "Heatmap_{source}_{rank}.pdf",
        ),
        parquet=pjoin(config["output_path"]["parquet"], "heatmap", "heatmap_{source}_{rank}.parquet"),
    benchmark:
        "benchmarks/{source}/heatmap_abundance_{rank}.tsv"
    wildcard_constraints:
        source="reads|contigs|kegg",
        rank="[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title_template=lambda w: get_text("TITLE_HEATMAP"),
        subtitle_template=lambda w: get_text("SUBTITLE_HEATMAP"),
        text_sample_all=lambda w: get_text("all"),
        top_n=config["plots"]["heatmap"]["top_n"],
        clust_method=config["plots"]["heatmap"]["clust_method"],
        distance_method=config["plots"]["heatmap"]["distance_method"],
        rank=lambda w: w.rank,
    script:
        os.path.abspath("workflow/scripts/plots/Heatmap_abundance.R")


rule plot_heatmap_deseq2:
    input:
        deseq_files=pjoin(
            config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"
        ),
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"
        ),
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "heatmap",
            "{source}",
            "Heatmap_deseq2_{source}_{rank}.pdf",
        ),
        parquet=pjoin(config["output_path"]["parquet"], "heatmap", "heatmap_deseq2_{source}_{rank}.parquet"),
    benchmark:
        "benchmarks/{source}/heatmap_deseq2_{rank}.tsv"
    wildcard_constraints:
        source="reads|contigs|kegg",
        rank="[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title_template=lambda w: get_text("TITLE_HEATMAP_DESEQ2"),
        subtitle_template=lambda w: get_text("SUBTITLE_HEATMAP_DESEQ2"),
        contrast=config["deseq2"]["contrast"],
        padj=config["plots"]["heatmap_deseq2"]["contrasts_values"]["padj_threshold"],
        lfc=config["plots"]["heatmap_deseq2"]["contrasts_values"]["lfc_threshold"],
        top_n=config["plots"]["heatmap_deseq2"]["top_n"],
        group_by=config["plots"]["heatmap_deseq2"]["group_by"],
        clust_method=config["plots"]["heatmap_deseq2"]["clust_method"],
        distance_method=config["plots"]["heatmap_deseq2"]["distance_method"],
        rank=lambda w: w.rank,
    script:
        os.path.abspath("workflow/scripts/plots/Heatmap_from_DESeq2.R")


rule plot_pca:
    input:
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"
        ),
        physico=config["input_path"]["physico_params"],
    output:
        pdf=pjoin(config["output_path"]["plots"], "pca", "{source}", "PCA_{source}_{rank}.pdf"),
        xlsx=pjoin(config["output_path"]["plots"], "pca", "pca_contributions_{source}_{rank}.xlsx"),
        parquet=pjoin(config["output_path"]["parquet"], "pca", "pca_{source}_{rank}.parquet"),
    benchmark:
        "benchmarks/{source}/pca_{rank}.tsv"
    wildcard_constraints:
        source="reads|contigs|kegg",
        rank="[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title_template=lambda w: get_text("TITLE_PCA"),
        subtitle_template=lambda w: get_text("SUBTITLE_PCA"),
        top_n=config["plots"]["pca"]["top_n"],
        dim_x=config["plots"]["pca"]["dim_x"],
        dim_y=config["plots"]["pca"]["dim_y"],
        physico_col=config["plots"]["pca"]["physico_col"],
        stand_col=lambda w: METRIC_PARAMS["normalized"][w.source]["stand_col"],
        point_size=config["plots"]["pca"]["point_size"],
        rank=lambda w: w.rank,
    script:
        os.path.abspath("workflow/scripts/plots/PCA.R")


rule plot_volcano:
    input:
        deseq_files=pjoin(
            config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"
        ),
        phyloseq_obj=pjoin(
            config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds"
        ),
    output:
        pdf=pjoin(
            config["output_path"]["plots"],
            "volcano",
            "{source}",
            "Volcano_{source}_{rank}.pdf",
        ),
        parquet=pjoin(
            config["output_path"]["parquet"],
            "volcano",
            "volcano_{source}_{rank}.parquet",
        ),
    benchmark:
        "benchmarks/{source}/volcano_{rank}.tsv"
    wildcard_constraints:
        source="reads|contigs|kegg",
        rank="[a-zA-Z0-9_]+",
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        title_template=lambda w: get_text("TITLE_VOLCANO"),
        subtitle_template=lambda w: get_text("SUBTITLE_VOLCANO"),
        contrast=config["deseq2"]["contrast"],
        padj=config["plots"]["volcano"]["contrasts_values"]["padj_threshold"],
        lfc=config["plots"]["volcano"]["contrasts_values"]["lfc_threshold"],
        top_n=config["plots"]["volcano"]["top_n"],
        rank=lambda w: w.rank,
    script:
        os.path.abspath("workflow/scripts/plots/Volcano.R")


rule plot_physico:
    input:
        physico_params=config["input_path"]["physico_params"]
    output:
        pdf=pjoin(config["output_path"]["plots"], "physico", "Physico_plot.pdf"),
        parquet=pjoin(
            config["output_path"]["parquet"], "physico", "physico_plot.parquet"
        )
    benchmark:
        "benchmarks/physico.tsv"
    conda:
        "envs/r_env.yaml"
    script:
        os.path.abspath("workflow/scripts/plots/Curves_physico_parameters.R")


# ==========================================================================
# DESEQ2 + PHYLOSEQ
# ==========================================================================

rule run_permanova:
    input:
        data=get_permanova_inputs,
        metadata=config["input_path"]["metadata"],
    output:
        pdf=pjoin(config["output_path"]["permanova"], "{source}", "Permanova_{source}.pdf"),
        xlsx=pjoin(config["output_path"]["permanova"], "{source}", "permanova_stat_{source}.xlsx"),
        parquet=pjoin(
            config["output_path"]["parquet"], "permanova", "permanova_{source}.parquet"
        )
    benchmark:
        "benchmarks/{source}/permanova.tsv"
    conda:
        "envs/r_env.yaml"
    params:
        shared=config["plots"]["shared"],
        stand_col = lambda wc: get_permanova_params(wc)["stand_col"],
        effect=config["permanova"]["effect"],
        distance_method=config["permanova"]["distance_method"],
        use_clr=config["permanova"]["use_clr"],
        permutation=config["permanova"]["permutation"],
    script:
        os.path.abspath("workflow/scripts/analysis/Permanova.R")


rule run_deseq2:
    input:
        data=lambda w: TREATMENT_SOURCES[f"{w.source}_counts"],
        metadata=config["input_path"]["metadata"],
    output:
        rds=pjoin(config["output_path"]["rds"], "{source}", "deseq2_{source}.rds"),
        parquet=pjoin(config["output_path"]["parquet"], "deseq2", "deseq2_{source}.parquet")
    benchmark:
        "benchmarks/{source}/deseq2.tsv"
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

rule run_phyloseq_other:
    input:
        data=lambda w: TREATMENT_SOURCES[f"{w.source}_normalized"],
        metadata=config["input_path"]["metadata"],
        taxonomy=lambda w: TAXONOMY[w.source],
    output:
        rds=pjoin(config["output_path"]["rds"], "{source}", "phyloseq_{source}_{rank}.rds")
    benchmark:
        "benchmarks/{source}/phyloseq_{source}_{rank}.tsv"
    wildcard_constraints:
        # Avoid a conflict (AmbiguousRuleException) with the KEGG rule
        source="reads|contigs"
    params:
        source=lambda w: w.source,
        stand_col=lambda w: METRIC_PARAMS["normalized"][w.source]["stand_col"],
        rank=lambda w: w.rank
    conda:
        "envs/r_env.yaml"
    script:
        os.path.abspath("workflow/scripts/analysis/Phyloseq.R")


rule run_phyloseq_kegg:
    input:
        data=lambda w: TREATMENT_SOURCES["kegg_normalized"],
        metadata=config["input_path"]["metadata"],
        taxonomy=lambda w: TAXONOMY["kegg"],
    output:
        rds=pjoin(config["output_path"]["rds"], "kegg", "phyloseq_kegg_{rank}.rds"),
    benchmark:
        "benchmarks/kegg/phyloseq_kegg_{rank}.tsv"
    wildcard_constraints:
        rank="level_3|gene_description" # Ensures {rank} only matches valid KEGG ranks
    params:
        source="kegg",
        stand_col=lambda w: METRIC_PARAMS["normalized"]["kegg"]["stand_col"],
        rank=lambda w: w.rank,  # Passes "ko" or "pathway" to R
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
        get_all_benchmarks
    output:
        csv="benchmarks/summary_benchmarks.csv"
    run:
        import pandas as pd
        from pathlib import Path

        records = []
        for file_path in input:
            path = Path(file_path)
            if path.exists():
                df = pd.read_csv(path, sep="\t")
                
                # Extract metadata from file structure
                df["rule"] = path.parent.name
                df["sample"] = path.stem.replace("_counted", "").replace("_filtered", "")
                records.append(df)

        if records:
            summary_df = pd.concat(records, ignore_index=True)
            
            # Reorder key benchmark columns if present
            primary_cols = ["rule", "sample", "s", "h:m:s", "max_rss", "cpu_time"]
            existing_primary = [c for c in primary_cols if c in summary_df.columns]
            other_cols = [c for c in summary_df.columns if c not in existing_primary]
            
            summary_df[existing_primary + other_cols].to_csv(output.csv, index=False)
        else:
            pd.DataFrame().to_csv(output.csv, index=False)
# ==========================================================================
# AUDIT
# ==========================================================================

ACTIVE_SOURCES = filter_active_sources(["reads", "contigs", "kegg"])
include: "audit/audit.smk"