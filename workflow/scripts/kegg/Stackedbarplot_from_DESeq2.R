################################################################################
# Project : "MicrobExplorer"
# Script: "Stackedbarplot from Kegg data DESeq2"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

library(phyloseq)
library(DESeq2)
library(ggplot2)
library(data.table)
library(viridis)

# ==========================================================================
# Configuration (Snakemake)
# ==========================================================================
res_path    <- snakemake[["input"]][["results_rds"]]  # "RDS_obj/deseq2_results.rds"
ps_path     <- snakemake[["input"]][["ps_kegg"]]      # "RDS_obj/ps_kegg_functional.rds"
pdf_out     <- snakemake[["output"]][["pdf"]]         # "Graphique/Barplot/Functional_Mosaics.pdf"
rds_out     <- snakemake[["output"]][["rds"]]         # "RDS_obj/processed_plot_data.rds"

# Thresholds from config
padj_thresh <- snakemake@config[["padj"]] %||% 0.05

# ==========================================================================
# 1. Data Preparation
# ==========================================================================
results_list <- readRDS(res_path)
ps_kegg      <- readRDS(ps_path)

# Relative abundance transformation
ps_rel <- transform_sample_counts(ps_kegg, function(x) x / sum(x))

# Optimized psmelt replacement using data.table
# This is much faster and memory-efficient for large KEGG tables
otu_dt  <- as.data.table(as(otu_table(ps_rel), "matrix"), keep.rownames = "OTU")
tax_dt  <- as.data.table(as(tax_table(ps_rel), "matrix"), keep.rownames = "OTU")
meta_dt <- as.data.table(as(sample_data(ps_rel), "data.frame"), keep.rownames = "Sample")

# Melt and merge
df_base <- melt(otu_dt, id.vars = "OTU", variable.name = "Sample", value.name = "Abundance")
df_base <- merge(df_base, tax_dt, by = "OTU")
df_base <- merge(df_base, meta_dt, by = "Sample")

# Cleanup memory
rm(otu_dt, tax_dt, ps_rel)
gc()

# ==========================================================================
# 2. Plotting Function
# ==========================================================================
generate_plot <- function(df, sig_ids, type_label, comp_name, ref_grp, test_grp, group_col) {
  
  # Filter data for the relevant groups
  df_relevant <- df[get(group_col) %in% c(ref_grp, test_grp)]
  
  # Identify Top 10 Significant KOs by abundance
  top_10 <- df_relevant[OTU %in% sig_ids, .(SumAb = sum(Abundance)), by = OTU][
    order(-SumAb)][1:10, OTU]
  
  # Labeling logic
  df_relevant[, Final_Category := "Non-Significatif"]
  df_relevant[OTU %in% sig_ids, Final_Category := "Autres KOs Significatifs"]
  
  # Add descriptions to Top 10
  # Note: assuming 'description' column exists in tax table
  df_relevant[OTU %in% top_10, Final_Category := paste0(OTU, ": ", substr(description, 1, 40), "...")]
  
  # Aggregate for plotting
  df_plot <- df_relevant[, .(Abundance = sum(Abundance)), by = .(Sample, Date, Final_Category)]
  
  # Colors
  unique_cats <- sort(setdiff(unique(df_plot$Final_Category), c("Autres KOs Significatifs", "Non-Significatif")))
  color_map <- c(
    setNames(viridis::turbo(length(unique_cats)), unique_cats),
    "Autres KOs Significatifs" = "grey60",
    "Non-Significatif"         = "grey10"
  )
  
  # Factor order
  df_plot[, Final_Category := factor(Final_Category, levels = c("Non-Significatif", "Autres KOs Significatifs", unique_cats))]
  
  # ggplot
  p <- ggplot(df_plot, aes(x = Sample, y = Abundance, fill = Final_Category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.9) +
    scale_fill_manual(values = color_map, name = "Function") +
    facet_grid(. ~ Date, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
      legend.text = element_text(size = 6),
      legend.position = "bottom",
      legend.direction = "vertical"
    ) +
    labs(
      title = paste(comp_name, "| Direction:", type_label),
      subtitle = paste("Comparison:", test_grp, "vs", ref_grp),
      y = "Relative Abundance"
    )
  
  return(p)
}

# ==========================================================================
# 3. Execution & Export
# ==========================================================================
pdf(pdf_out, width = 12, height = 9)

for (nom_comp in names(results_list)) {
  res <- as.data.table(results_list[[nom_comp]])
  
  # Detection of grouping type (Date vs Digesteur)
  is_date <- "Date_Test" %in% names(res)
  t_act   <- if(is_date) res$Date_Test[1] else res$Test_Group[1]
  t_ref   <- if(is_date) res$Date_Ref[1]  else res$Ref_Group[1]
  grp_col <- if(is_date) "Date_Raw"       else "Digesteur"
  
  # Filter sig KOs based on config padj
  ids_over  <- res[padj <= padj_thresh & log2FoldChange > 0, KO_Number]
  ids_under <- res[padj <= padj_thresh & log2FoldChange < 0, KO_Number]
  
  if(length(ids_over) > 0) {
    print(generate_plot(df_base, ids_over, "Over-represented", nom_comp, t_ref, t_act, grp_col))
  }
  if(length(ids_under) > 0) {
    print(generate_plot(df_base, ids_under, "Under-represented", nom_comp, t_ref, t_act, grp_col))
  }
}

dev.off()

# Save processed data table for future use in Shiny
saveRDS(df_base, rds_out)

message("✓ PDF generated: ", pdf_out)
message("✓ Plotting data saved: ", rds_out)