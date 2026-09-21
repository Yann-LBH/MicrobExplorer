# ==============================================================================
# PROJECT   : MicrobExplorer
# SCRIPT    : utils_stackedbarplot_from_DESeq2.R
# PURPOSE   : A utility that manages the operation of stackedbarplot_from_DESeq2
# AUTHOR    : Yann Le Bihan
# DATE      : 2026-09-03
# LINK      : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

# Libraries CRAN
library(data.table)
library(ggplot2)
library(viridis)

# Libraries Bioconductor
library(phyloseq)

# ==========================================================================
# 1. Charger les résultats DESeq2
# ==========================================================================
load_deseq2_results <- function(deseq_paths) {
  #' Load DESeq2 results from the new Master RDS structure
  #'
  #' @param deseq_paths Named list or named vector with file paths
  #'        Names should be contrast types (e.g., "ref", "date", "combo")
  #' @return List of data.tables with DESeq2 results

  results <- list()

  if (!is.list(deseq_paths)) {
    deseq_paths <- as.list(deseq_paths)
  }

  for (contrast_type in names(deseq_paths)) {
    file_path <- deseq_paths[[contrast_type]]

    if (is.na(file_path) || is.null(file_path)) {
      warning(sprintf("Skipping %s: path is NULL or NA", contrast_type))
      next
    }

    if (file.exists(file_path)) {
      tryCatch(
        {
          res_master <- readRDS(file_path)

          # 🟢 Adaptation à la nouvelle structure : On cible directement le slot $dt
          if (is.list(res_master) && !is.null(res_master$dt)) {
            # Si le type demandé (ex: 'ref') existe dans les tables pré-calculées
            if (contrast_type %in% names(res_master$dt)) {
              res_dt <- as.data.table(res_master$dt[[contrast_type]])
            } else {
              # Fallback : Si l'utilisateur passe directement le sous-slot ou une table globale
              res_dt <- as.data.table(res_master$dt)
            }
          } else if (is.data.frame(res_master) || inherits(res_master, "data.table")) {
            # Rétrocompatibilité au cas où un tableau brut soit passé
            res_dt <- as.data.table(res_master)
          } else {
            stop(sprintf("Unexpected Master RDS structure for %s", contrast_type))
          }

          # Standardisation uniforme de la colonne d'identifiant pour la suite de l'utilitaire
          if (!"Feature_ID" %in% names(res_dt) && "KO_Number" %in% names(res_dt)) {
            setnames(res_dt, "KO_Number", "Feature_ID")
          } else if (!"Feature_ID" %in% names(res_dt)) {
            setnames(res_dt, 1, "Feature_ID")
          }

          results[[contrast_type]] <- res_dt
          message(sprintf("✓ Loaded %s contrast table: %d features", contrast_type, nrow(res_dt)))
        },
        error = function(e) {
          warning(sprintf("Failed to load %s: %s", contrast_type, e$message))
        }
      )
    } else {
      warning(sprintf("File not found for %s: %s", contrast_type, file_path))
    }
  }

  if (length(results) == 0) {
    stop("No DESeq2 results could be loaded")
  }

  return(results)
}

# ==========================================================================
# 2. Transformer phyloseq en data.table pour relative abundance
# ==========================================================================
phyloseq_to_dt <- function(ps_object) {
  #' Convert phyloseq object to data.table format
  #'
  #' @param ps_object phyloseq object (should already be relative abundance)
  #' @return data.table with columns: OTU, Sample, Abundance, Metadata...

  otu_dt <- as.data.table(as(otu_table(ps_object), "matrix"), keep.rownames = "OTU")
  tax_dt <- as.data.table(as(tax_table(ps_object), "matrix"), keep.rownames = "OTU")
  meta_dt <- as.data.table(as(sample_data(ps_object), "data.frame"), keep.rownames = "Sample")

  df <- melt(otu_dt, id.vars = "OTU", variable.name = "Sample", value.name = "Abundance")
  df <- merge(df, tax_dt, by = "OTU")
  df <- merge(df, meta_dt, by = "Sample")

  return(df)
}

# ==========================================================================
# 3. Obtenir les features significatives
# ==========================================================================
get_significant_features <- function(DESEQ_RESULTS, PADJ_THRESHOLD, LFC_THRESHOLD, direction) {
  #' Get significant features from DESeq2 results
  #'
  #' @param DESEQ_RESULTS data.table with DESeq2 results
  #' @param PADJ_THRESHOLD p-adjusted threshold (default 0.05)
  #' @param LFC_THRESHOLD log2 fold change threshold (default 1.3)
  #' @param direction "up" or "down" for over/under-represented
  #' @return Vector of feature IDs

  # Sécurité : mapper le nom de la colonne d'identifiant si Feature_ID s'appelle KO_Number
  id_col <- if ("Feature_ID" %in% names(DESEQ_RESULTS)) "Feature_ID" else "KO_Number"

  if (direction == "up") {
    features <- DESEQ_RESULTS[padj <= PADJ_THRESHOLD & log2FoldChange >= LFC_THRESHOLD][[id_col]]
  } else if (direction == "down") {
    features <- DESEQ_RESULTS[padj <= PADJ_THRESHOLD & log2FoldChange <= -LFC_THRESHOLD][[id_col]]
  } else {
    stop("direction must be 'up' or 'down'")
  }

  return(features)
}

# ==========================================================================
# 4. Générer la palette de couleurs
# ==========================================================================
generate_color_palette <- function(top_features, palette_name = "turbo") {
  #' @param top_features Vector of top feature names (Display_Names)
  #' @param palette_name Character string matching a viridis function (e.g., "turbo", "magma", "plasma")
  #' @return Named vector of colors
  
  # Sécurité : Si la palette demandée n'existe pas dans viridis, on se rabat sur turbo
  if (!palette_name %in% c("turbo", "viridis", "magma", "plasma", "inferno", "cividis")) {
    palette_name <- "turbo"
  }

  # Sécurisation : extraction des noms uniques pour éviter les doublons
  top_features_unique <- unique(top_features)

  # Appel dynamique de la fonction du package viridis (ex: viridis::turbo(n))
  colors <- do.call(get(palette_name, envir = asNamespace("viridis")), list(length(top_features)))
  names(colors) <- top_features_unique
  
  return(colors)
}

# ==========================================================================
# 5. Préparer les données pour le plot (Version Corrigée et Sécurisée)
# ==========================================================================
prepare_plot_data <- function(df_abundance, sig_features, TOP_N, group_column = "group", rank_column = NULL, comp_name = NULL, time_column = "date") {
  #' Prepare data for stacked barplot with auto-filtering based on contrast elements

  df_local <- copy(df_abundance)

  # 🔍 DÉTECTION ET FILTRAGE DYNAMIQUE SELON LE CONTRASTE
  if (!is.null(comp_name) && grepl("_vs_", comp_name)) {
    elements_in_contrast <- unlist(strsplit(comp_name, "_vs_"))
    elements_in_contrast <- trimws(elements_in_contrast)
    
    if (time_column %in% names(df_local)) {
      df_local[, time_var_str := as.character(get(time_column))]
    } else {
      df_local[, time_var_str := character(0)]
    }
    
    time_values_str <- unique(df_local$time_var_str)
    
    # Mode Temporel
    if (length(time_values_str) > 0 && any(elements_in_contrast %in% time_values_str)) {
      df_local <- df_local[time_var_str %in% elements_in_contrast]
      df_local[, Group_Plot_Var := paste0(get(group_column), "_", time_var_str)]
      target_group_col <- "Group_Plot_Var"
      message(sprintf("    [Plot Prep] Automated TIME mode active for %s. Samples retained: %d", comp_name, nrow(df_local)))
    } else {
      # Mode Condition/Ref/Combo
      if (group_column %in% names(df_local)) {
        group_values_str <- as.character(unique(df_local[[group_column]]))
        if (all(elements_in_contrast %in% group_values_str)) {
          df_local <- df_local[as.character(get(group_column)) %in% elements_in_contrast]
        }
      }
      target_group_col <- group_column
      message(sprintf("    [Plot Prep] Automated CONDITION mode active for %s. Samples retained: %d", comp_name, nrow(df_local)))
    }
  } else {
    target_group_col <- group_column
  }

  if ("time_var_str" %in% names(df_local)) df_local[, time_var_str := NULL]

  if (nrow(df_local) == 0) {
    warning(sprintf("Warning: Zero rows remaining for comparison '%s' after dynamic filtering.", comp_name))
    return(data.table(Group_Var = character(0), Category = factor(character(0)), Abundance = numeric(0)))
  }

  match_rank_col <- names(df_local)[tolower(names(df_local)) == tolower(rank_column)]

  if (length(match_rank_col) > 0 && match_rank_col[1] != "OTU") {
    # 1. Cas Taxonomie classique (ex: Species, Genus, Family...) pour contigs / reads
    val_rank <- as.character(df_local[[match_rank_col[1]]])
    
    # Nettoyage des chaînes vides ou NA
    val_rank[is.na(val_rank) | val_rank == "" | val_rank == "NA"] <- "Unassigned"
    
    df_local[, Display_Name := val_rank]

  } else if ("combined_label" %in% names(df_local) && any(!is.na(df_local$combined_label))) {
    # 2. Cas KEGG (KO | description)
    df_local[, Display_Name := as.character(combined_label)]

  } else {
    # 3. Fallback sur l'OTU / ID du contig
    df_local[, Display_Name := as.character(OTU)]
  }

  # 🟢 NETTOYAGE & TRONCATURE DU LIBELLÉ
  df_local[is.na(Display_Name) | Display_Name == "" | Display_Name == "NA", Display_Name := "Unassigned"]
  
  df_local[!grepl("unassigned", Display_Name, ignore.case = TRUE), 
           Display_Name := ifelse(nchar(Display_Name) > 45, paste0(substr(Display_Name, 1, 42), "..."), Display_Name)]

  df_sig <- df_local[OTU %in% sig_features]

  top_features <- character(0)
  if (nrow(df_sig) > 0) {
    top_dt <- df_sig[, .(Total_Abundance = sum(Abundance)), by = Display_Name][order(-Total_Abundance)]
    top_features <- head(top_dt$Display_Name, TOP_N)
  }

  df_local[, Category := "Not significant"]
  df_local[OTU %in% sig_features, Category := "Significant (Other)"]
  df_local[Display_Name %in% top_features & OTU %in% sig_features, Category := Display_Name]

  df_plot <- df_local[, .(Abundance = sum(Abundance)),
    by = .(Group_Var = get(target_group_col), Category)
  ]
  
  df_plot[, Total_Group_Abund := sum(Abundance), by = Group_Var]
  df_plot[Total_Group_Abund > 0, Abundance := Abundance / Total_Group_Abund]

  factor_levels <- c("Not significant", "Significant (Other)", sort(top_features))
  df_plot[, Category := factor(Category, levels = factor_levels)]

  return(df_plot)
}

# ==========================================================================
# 6. Créer le plot stackedbarplot (Version Corrigée et Sécurisée)
# ==========================================================================
create_stackedbarplot <- function(df_plot, title, subtitle, feature_colors, group_label = "condition") {
  #' @param df_plot data.table processed by prepare_plot_data
  #' @param title Plot title
  #' @param subtitle Plot subtitle
  #' @param feature_colors Named vector generated by generate_color_palette
  #' @param group_label Axis label name
  #' @return ggplot object

  color_map <- c(
    feature_colors,
    "Significant (Other)" = "grey50",
    "Not significant"     = "#000000"
  )

  p <- ggplot(df_plot, aes(x = Group_Var, y = Abundance, fill = Category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.75) +
    scale_fill_manual(values = color_map, name = "Features") +
    scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "right", 
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9)
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = group_label,
      y = "Mean Relative Abundance"
    )

  return(p)
}

# ==========================================================================
# Export functions
# ==========================================================================
export_results <- function(df_plot, output_path) {
  tryCatch(
    {
      arrow::write_parquet(df_plot, output_path)
      message(sprintf("✓ Results exported to: %s", output_path))
    },
    error = function(e) {
      warning(sprintf("Failed to export results: %s", e$message))
    }
  )
}