rule project_audit:
    input:
        metadata="data/metadata/metadata.xlsx",
        reads_data=expand(READS_TREATMENT + "4.Annotated/annotated_{sample}_reads.tsv", sample=SAMPLES),
        contigs_data=expand(CONTIGS_TREATMENT + "6.Annotated/annotated_{sample}_contigs.tsv", sample=SAMPLES),
        kegg_data=expand(KEGG_TREATMENT + "023.Annotated/annotated_{sample}_kegg.tsv", sample=SAMPLES),
        ps_kegg="results/rds/kegg/phyloseq_kegg.rds"
    output:
        report=config["output_path"]["audit"]
    log:
        "logs/audit/project_audit.log"
    conda:
        "../envs/r_env.yaml"
    script:
        "audit.R"