library(readxl)
library(data.table)
library(rmarkdown)
library(phyloseq)

# Fetch Snakemake variables
meta_path   <- snakemake@input[["metadata"]]
reads_files <- snakemake@input[["reads_data"]]
contigs_files <- snakemake@input[["contigs_data"]]
kegg_files  <- snakemake@input[["kegg_data"]]
ps_kegg_path <- snakemake@input[["ps_kegg"]]
output_html <- snakemake@output[["report"]]

# ==========================================================================
# AUDIT STEP 1: Sample Completeness and Tracking
# ==========================================================================
meta <- as.data.table(read_excel(meta_path))
expected_samples <- meta$sample_id

# Helper function to extract sample IDs from file paths
get_samples_from_paths <- function(files) {
  sapply(files, function(f) {
    # Match any sample ID from metadata present in the filename
    matched <- expected_samples[sapply(expected_samples, function(sid) grepl(sid, basename(f)))]
    if(length(matched) == 0) return(NA) else return(matched[1])
  })
}

reads_samples   <- na.omit(get_samples_from_paths(reads_files))
contigs_samples <- na.omit(get_samples_from_paths(contigs_files))
kegg_samples    <- na.omit(get_samples_from_paths(kegg_files))

# Check if any sample went missing during the processing steps
missing_reads   <- setdiff(expected_samples, reads_samples)
missing_contigs <- setdiff(expected_samples, contigs_samples)
missing_kegg    <- setdiff(expected_samples, kegg_samples)

# ==========================================================================
# AUDIT STEP 2: Multi-omic Correlation Check
# ==========================================================================
# We check if features/counts trends are correlated across datatypes per sample
audit_summary <- data.table(sample_id = expected_samples)

# Get total features or sums per sample type to build a trend matrix
audit_summary[, Total_Reads_Features := sapply(reads_files, function(f) nrow(fread(f, select=1)))]
audit_summary[, Total_Contigs_Features := sapply(contigs_files, function(f) nrow(fread(f, select=1)))]
audit_summary[, Total_KEGG_Features := sapply(kegg_files, function(f) nrow(fread(f, select=1)))]

# ==========================================================================
# AUDIT STEP 3: Phyloseq Object Validation
# ==========================================================================
ps <- readRDS(ps_kegg_path)
ps_samples <- sample_names(ps)
ps_taxa_count <- ntaxa(ps)
is_ps_corrupted <- any(duplicated(taxa_names(ps)))

# ==========================================================================
# GENERATE HTML REPORT
# ==========================================================================

# 1. Pre-calculate the Sample Tracking status HTML to avoid inline syntax syntax errors
sample_status_html <- if (length(missing_reads) == 0 && length(missing_contigs) == 0 && length(missing_kegg) == 0) {
  paste0("<div class='card success'><strong>PASSED:</strong> All expected samples (", length(expected_samples), ") are present across Reads, Contigs, and KEGG datasets.</div>")
} else {
  paste0("<div class='card alert'><strong>WARNING:</strong> Missing samples detected!<br>",
         "Missing in Reads: ", paste(missing_reads, collapse=", "), "<br>",
         "Missing in Contigs: ", paste(missing_contigs, collapse=", "), "<br>",
         "Missing in KEGG: ", paste(missing_kegg, collapse=", "), "</div>")
}

# Create a dynamic HTML report summarizing findings and potential alerts
html_content <- paste0("
<html>
<head>
    <style>
        body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; color: #333; }
        h1 { color: #2c3e50; border-bottom: 2px solid #ecf0f1; padding-bottom: 10px; }
        h2 { color: #16a085; margin-top: 30px; }
        .card { background: #f8f9fa; border-left: 5px solid #3498db; padding: 15px; margin: 10px 0; border-radius: 4px; }
        .alert { border-left-color: #e74c3c; background: #fdf2f2; }
        .success { border-left-color: #2ecc71; background: #f4fbf7; }
        table { width: 100%; border-collapse: collapse; margin-top: 15px; }
        th, td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        th { background-color: #2c3e50; color: white; }
    </style>
</head>
<body>
    <h1>MicrobExplorer - Pipeline Global Audit Report</h1>
    <p>Generated on: ", Sys.time(), "</p>
    
    <h2>1. Sample Tracking & Integrity</h2>
    ", sample_status_html, "

    <h2>2. Phyloseq Sanity Check</h2>
    <div class='card ", if(!is_ps_corrupted) "success" else "alert", "'>
        <strong>Phyloseq Object status:</strong> Total features integrated = ", ps_taxa_count, "<br>",
        "<strong>Duplicated Taxa Names Check:</strong> ", if(!is_ps_corrupted) "CLEAN (0 duplicates found)" else "CRITICAL ERROR: Duplicates remaining!", "<br>",
        "<strong>Samples successfully loaded in RDS:</strong> ", length(ps_samples), " / ", length(expected_samples), "
    </div>

    <h2>3. Multi-omic Data Density Matrix</h2>
    <p>Comparison of feature counts extracted per sample across all omic layers.</p>
    <table>
        <tr><th>Sample ID</th><th>Reads Rows</th><th>Contigs Rows</th><th>KEGG Rows</th></tr>",
        paste(sapply(1:nrow(audit_summary), function(i) {
            paste0("<tr><td>", audit_summary$sample_id[i], "</td><td>", audit_summary$Total_Reads_Features[i], "</td><td>", audit_summary$Total_Contigs_Features[i], "</td><td>", audit_summary$Total_KEGG_Features[i], "</td></tr>")
        }), collapse=""), "
    </table>
</body>
</html>
")

writeLines(html_content, output_html)