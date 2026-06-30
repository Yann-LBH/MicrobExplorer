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

# Inputs / Outputs
DATA     <- as.character(snakemake@input[["data"]])
METADATA <- as.character(snakemake@input[["metadata"]])
RDS      <- as.character(snakemake@output[["rds"]])

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
# 2. Détection Automatique du Type de Données (KEGG vs Taxonomie)
# ==========================================================================

 if ("cpm" %in% names(all_data)) {
  # --- Configuration MODE TAXONOMIE (Reads) ---
  # ✅ FIXED: Use read_id as row identifier for the reads dataset
  id_col   <- "read_id" 
  tax_cols <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_cols <- intersect(tax_cols, names(all_data))
  message("🦠 Mode detected : Taxonomic (Reads)")
}  else if ("rpkm" %in% names(all_data)) {
  # --- Configuration MODE TAXONOMIE (Contigs) ---
  # ✅ FIXED: Use contig_id as row identifier, matching the lowercase python step 6
  id_col   <- "contig_id"
  tax_cols <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  tax_cols <- intersect(tax_cols, names(all_data))
  message("🦠 Mode detected : Taxonomic (Contigs)")
} else {
  # --- Configuration MODE KEGG ---
  id_col   <- "ko"
  tax_cols <- intersect(c("ec_number", "level_1", "level_2", "level_3", "gene_description"), names(all_data))
  message("🧬 Mode detected : Functional (KEGG)")

}

# Sécurité critique : On vérifie que la colonne d'abondance demandée existe
if (!VALUE_COL %in% names(all_data)) {
  stop(sprintf("Error : The abundance column [%s] does not exist in these files.", VALUE_COL))
}

# ==========================================================================
# 3. Construction des Composants Phyloseq
# ==========================================================================

# --- A. OTU TABLE (Matrice d'abondance) ---
formula_str <- as.formula(paste(id_col, "~ sample_id"))
otu_dt <- dcast(all_data, formula_str, value.var = VALUE_COL, fun.aggregate = sum, fill = 0)

otu_mat <- as.matrix(otu_dt, rownames = id_col)
mode(otu_mat) <- "numeric"

# --- B. TAX TABLE (Table de classification) ---
tax_dt  <- unique(all_data[, c(id_col, tax_cols), with = FALSE])
tax_mat <- as.matrix(tax_dt, rownames = id_col)

# --- C. SAMPLE DATA (Métadonnées) ---
sample_df <- as.data.frame(meta_dt[sample_id %in% colnames(otu_mat)])
rownames(sample_df) <- sample_df$sample_id

# ==========================================================================
# 4. Assemblage, Normalisation et Sauvegarde
# ==========================================================================
ps_final <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(sample_df)
)

# Normalisation automatique si on part de reads bruts (évite les biais de séquençage)
if (VALUE_COL == "reads_count") {
  ps_final <- transform_sample_counts(ps_final, function(x) x / sum(x))
  message("✓ Normalisation effectuée : transformation des counts bruts en abondances relatives.")
}

saveRDS(ps_final, RDS)
message("✓ Objet Phyloseq RDS créé avec succès : ", RDS)
