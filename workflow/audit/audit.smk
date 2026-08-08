rule project_audit:
    input:
        metadata="data/metadata/metadata.xlsx",
        reads_data=expand(READS_TREATMENT + "3.CPM/cpm_{sample}_reads.tsv", sample=SAMPLES) if "reads" in ACTIVE_SOURCES else [],
        contigs_data=expand(CONTIGS_TREATMENT + "5.Union/union_{sample}_contigs.tsv", sample=SAMPLES) if "contigs" in ACTIVE_SOURCES else [],
        kegg_data=expand(KEGG_TREATMENT + "3.RPKM/rpkm_{sample}_kegg.tsv", sample=SAMPLES) if "kegg" in ACTIVE_SOURCES else [],
        ps_kegg="results/rds/kegg/phyloseq_kegg.rds" if "kegg" in ACTIVE_SOURCES else []
    output:
        report=config["output_path"]["audit"]
    log:
        "logs/audit/project_audit.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "audit.R"