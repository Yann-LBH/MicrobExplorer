library(data.table)
library(phyloseq)
library(arrow)
library(readxl)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================
data <- snakemake[["input"]][["data_dir"]]
metadata <- snakemake[["input"]][["metadata"]]
rds_save <- snakemake[["output"]][["rds"]]

# ==========================================================================
# 1. Data Import
# ==========================================================================
# Metadata
meta_dt <- as.data.table(read_xlsx(metadata))
setkey(meta_dt, sample_id)

# KO Data - Fast multi-file loading
files <- list.files(data, pattern = "annotated_.*\\.tsv", full.names = TRUE)

df_list <- lapply(files, function(f) {
  # Fast ID extraction
  sample_nm <- gsub("annotated_agreg_|\\.tsv", "", basename(f))
  dt <- fread(f, select = c(
    "ko", "standardization", "ec_number",
    "level_1", "level_2", "level_3", "gene_description"
  ))
  dt[, sample_id := sample_nm]
  return(dt)
})

all_kos <- rbindlist(df_list)

# ==========================================================================
# 2. Matrix Construction
# ==========================================================================

# --- A. OTU TABLE (Counts per KO) ---
# Data.table dcast is significantly faster than pivot_wider
otu_dt <- dcast(all_kos, ko ~ sample_id,
  value.var = "standardization",
  fun.aggregate = sum,
  fill = 0
)

otu_mat <- as.matrix(otu_dt, rownames = "ko")

# --- B. TAX TABLE (KEGG Hierarchy) ---
tax_dt <- unique(all_kos[, .(ko, ec_number, level_1, level_2, level_3, gene_description)])
tax_mat <- as.matrix(tax_dt, rownames = "ko")

# --- C. Metadata preparation ---
sample_df <- as.data.frame(meta_dt)
rownames(sample_df) <- sample_df$sample_id

# ==========================================================================
# 3. Phyloseq Initialization & Export
# ==========================================================================
ps_kegg <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

# --- Export to RDS ---
saveRDS(ps_kegg, rds_save)

message("✓ Phyloseq RDS object created.")
