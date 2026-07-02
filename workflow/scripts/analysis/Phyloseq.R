################################################################################
# Project : "MicrobExplorer"
# Script  : "Analysis : Phyloseq"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================
library(data.table)
library(phyloseq)
library(readxl)

# Inputs
DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])[1]

# Outputs
RDS      <- as.character(snakemake@output[["rds"]])[1]

# Paramètres (avec valeurs par défaut au cas où)
VALUE_COL <- tolower(as.character(snakemake@params[["value_col"]]))[1]

# ==========================================================================
# 1. Chargement des métadonnées et des fichiers TSV
# ==========================================================================
meta_dt <- as.data.table(read_xlsx(METADATA))
meta_dt[, sample_id := as.character(sample_id)]
setkey(meta_dt, sample_id)

df_list <- lapply(DATA, function(f) {
  file_name <- basename(f)
  matched_sample <- meta_dt[sapply(sample_id, function(sid) grepl(sid, file_name)), sample_id]
  
  if (length(matched_sample) == 0 || is.na(matched_sample)) return(NULL)
  
  dt <- fread(f, showProgress = FALSE)
  if (nrow(dt) == 0) return(NULL)
  
  dt[, sample_id := matched_sample]
  return(dt)
})

all_data <- rbindlist(Filter(Negate(is.null), df_list), use.names = TRUE, fill = TRUE)

# ==========================================================================
# 2. Automatic Data Type Detection (Reads vs Contigs vs KEGG)
# ==========================================================================
if ("cpm" %in% names(all_data)) {
  id_col   <- "read_id" 
  tax_cols <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_cols <- intersect(tax_cols, names(all_data))
  message("🦠 Mode detected: Taxonomic (Reads)")
} else if ("rpkm" %in% names(all_data)) {
  id_col   <- "contig_id"
  tax_cols <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_cols <- intersect(tax_cols, names(all_data))
  message("🦠 Mode detected: Taxonomic (Contigs)")
} else {
  id_col   <- "ko"
  tax_cols <- intersect(c("pathway", "description", "ec_number", "level_1", "level_2", "level_3", "gene_description"), names(all_data))
  message("🧬 Mode detected: Functional (KEGG)")
}

# Security check for mandatory columns
if (is.na(id_col) || !id_col %in% names(all_data)) {
  stop("Error: The dynamic ID column could not be resolved or is missing from data.")
}
if (!VALUE_COL %in% names(all_data)) {
  stop(sprintf("Error: The abundance column [%s] does not exist in these files.", VALUE_COL))
}

# Format IDs as character and clean missing values
all_data[, (id_col) := as.character(get(id_col))]
all_data <- all_data[!is.na(get(id_col)) & get(id_col) != "" & get(id_col) != "NA"]

# ==========================================================================
# 3. Build Phyloseq Components
# ==========================================================================

# --- A. OTU TABLE ---
# fun.aggregate = sum handles multiple identical KOs per sample perfectly
formula_str <- as.formula(paste(id_col, "~ sample_id"))
otu_dt <- dcast(all_data, formula_str, value.var = VALUE_COL, fun.aggregate = sum, fill = 0)

otu_mat <- as.matrix(otu_dt, rownames = id_col)
mode(otu_mat) <- "numeric"

# --- B. TAX TABLE (FIXED FOR KEGG DUPLICATES) ---
# Select ID and tax columns, then use unique(..., by = id_col) to keep exactly ONE row per ID
tax_dt  <- unique(all_data[, c(id_col, tax_cols), with = FALSE], by = id_col)
tax_mat <- as.matrix(tax_dt, rownames = id_col)

# --- C. SAMPLE DATA ---
sample_df <- as.data.frame(meta_dt[sample_id %in% colnames(otu_mat)])
rownames(sample_df) <- sample_df$sample_id

# ==========================================================================
# 4. Assemble and Save Object
# ==========================================================================
ps_final <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

saveRDS(ps_final, RDS)
message("✓ Phyloseq RDS object successfully created: ", RDS)