# ==============================================================================
# PROJECT   : MicrobExplorer
# SCRIPT    : utils_io.R
# PURPOSE   : Centralize data, metadata, taxonomy and count matrix loading
# AUTHOR    : Yann Le Bihan
# DATE      : 2026-09-03
# LINK      : https://github.com/Yann-LBH/MicrobExplorer
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(tools)
  library(readxl)
})

# Colonnes obligatoires dans tout fichier de metadonnees consomme par le pipeline
REQUIRED_METADATA_COLS <- c("sample_id", "name", "date")
 
# ==========================================================================
# Helpers internes
# ==========================================================================
 
#' Escape a string so it can be safely embedded in a regular expression.
#' @param x Character vector to escape.
#' @return Character vector with regex-special characters escaped.
.regex_escape <- function(x) {
  gsub("([.\\+*?^$()\\[\\]{}|\\\\])", "\\\\\\1", x, perl = TRUE)
}
 
#' Find, among valid_ids, which ones appear in file_name as a whole token
#' (bounded by start/end of string or a non-alphanumeric separator).
#' Returns the longest match (most specific) when several ids match, and
#' raises a warning if several ids of maximal length match ambiguously.
#' @param file_name Character scalar, the file name to search into.
#' @param valid_ids Character vector of candidate sample_id values.
#' @return A single matched sample_id, or NA_character_ if none matched.

.match_sample_id <- function(file_name, valid_ids) {
    file_base <- tools::file_path_sans_ext(file_name)
    matches <- valid_ids[vapply(valid_ids, function(id) {
        pattern <- paste0("(^|[^A-Za-z0-9])", .regex_escape(id), "([^A-Za-z0-9]|$)")
        grepl(pattern, file_base)
    }, logical(1))]
    
    if (length(matches) == 0) return(NA_character_)
    
    max_len <- max(nchar(matches))
    best_matches <- matches[nchar(matches) == max_len]
    
    if (length(best_matches) > 1) {
        warning(sprintf(
        "Ambiguous sample_id match for file '%s': candidates [%s] have equal specificity. Using '%s'.",
        file_name, paste(best_matches, collapse = ", "), best_matches[1]
        ))
    }
    
    best_matches[1]
}
 
# ==========================================================================
# Chargement des metadonnees
# ==========================================================================
 
#' Loads a metadata file (Excel or delimited text) and ensures it satisfies
#' the minimal contract required by the pipeline: presence of sample_id,
#' name and date columns, uniqueness of sample_id, and character typing of
#' the identifying columns.
#'
#' @param path Path to the metadata file (.xlsx, .xls, .csv, .tsv, ...).
#' @return data.table keyed on sample_id, with sample_id/name/date as character.

load_metadata <- function(path) {
    if (!file.exists(path)) {
        stop(sprintf("Metadata file not found: %s", path))
    }
    
    ext <- tolower(tools::file_ext(path))
    dt <- if (ext %in% c("xlsx", "xls")) {
        data.table::as.data.table(readxl::read_excel(path))
    } else {
        data.table::fread(path, showProgress = FALSE)
    }
    
    missing_cols <- setdiff(REQUIRED_METADATA_COLS, names(dt))
    if (length(missing_cols) > 0) {
        stop(sprintf(
        "Metadata file '%s' is missing mandatory column(s): %s",
        basename(path), paste(missing_cols, collapse = ", ")
        ))
    }
    
    # Conversion securisee des colonnes identifiantes en character
    cols_to_convert <- intersect(c("sample_id", "name", "date"), names(dt))
    dt[, (cols_to_convert) := lapply(.SD, as.character), .SDcols = cols_to_convert]
    
    # sample_id ne doit contenir ni NA/vide, ni doublon
    empty_ids <- is.na(dt$sample_id) | dt$sample_id == ""
    if (any(empty_ids)) {
        stop(sprintf(
        "Metadata file '%s' contains %d row(s) with a missing/empty 'sample_id'.",
        basename(path), sum(empty_ids)
        ))
    }
    
    dupes <- unique(dt$sample_id[duplicated(dt$sample_id)])
    if (length(dupes) > 0) {
        stop(sprintf(
        "Metadata file '%s' contains duplicated sample_id(s): %s",
        basename(path), paste(dupes, collapse = ", ")
        ))
    }
    
    # name/date: on tolere les NA (pas bloquant pour le pipeline) mais on avertit
    empty_name <- is.na(dt$name) | dt$name == ""
    if (any(empty_name)) {
        warning(sprintf(
        "Metadata file '%s' has %d row(s) with a missing/empty 'name'.",
        basename(path), sum(empty_name)
        ))
    }
    empty_date <- is.na(dt$date) | dt$date == ""
    if (any(empty_date)) {
        warning(sprintf(
        "Metadata file '%s' has %d row(s) with a missing/empty 'date'.",
        basename(path), sum(empty_date)
        ))
    }
    
    data.table::setkey(dt, sample_id)
    dt
}
 
# ==========================================================================
# Chargement des TSV et association au sample_id
# ==========================================================================
 
#' Dynamically loads a list of TSV files and associates each one with a
#' sample_id from the metadata, based on the file name. Matching is done on
#' whole tokens (bounded by non-alphanumeric separators or string
#' boundaries) to avoid partial/substring false positives (e.g. "D1"
#' matching inside "STD182"). When several sample_id match a file name with
#' equal specificity, the ambiguity is logged via warning() and the first
#' candidate is used.
#'
#' @param paths Vector of paths to the TSV files.
#' @param meta_dt data.table of metadata containing the sample_id column
#'   (typically the output of load_metadata()).
#' @param select_cols Optional vector of columns to read (significantly
#'   speeds up reading on large files). 'sample_id' is automatically
#'   excluded from the read since it is assigned after matching.
#' @return Consolidated data.table with a `sample_id` column.
load_tsv_dir_dynamic <- function(paths, meta_dt, select_cols = NULL) {
 
    valid_ids <- as.character(meta_dt$sample_id)
    
    rows <- lapply(paths, function(f) {
        if (!file.exists(f)) {
        message(sprintf("File not found, skipped: %s", f))
        return(NULL)
        }
    
        file_name <- basename(f)
        matched_sample <- .match_sample_id(file_name, valid_ids)
    
        if (is.na(matched_sample)) {
        message(sprintf("File ignored (no sample_id found in the filename): %s", file_name))
        return(NULL)
        }
    
        cols_to_read <- if (!is.null(select_cols)) setdiff(select_cols, "sample_id") else NULL
    
        dt <- tryCatch({
        if (is.null(cols_to_read)) {
            data.table::fread(f, showProgress = FALSE)
        } else {
            data.table::fread(f, select = cols_to_read, showProgress = FALSE)
        }
        }, error = function(e) {
        cols_str <- if (!is.null(cols_to_read)) paste(cols_to_read, collapse = ", ") else "ALL"
        stop(sprintf(
            "Error reading file '%s': unable to load columns [%s].\nOriginal message: %s",
            file_name, cols_str, e$message
        ))
        })
    
        if (nrow(dt) == 0) return(NULL)
    
        dt[, sample_id := matched_sample]
        dt
    })
    
    dt_final <- data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
    
    if (nrow(dt_final) == 0 || !"sample_id" %in% names(dt_final)) {
        stop("CRITICAL ERROR: No TSV file could be associated with a sample_id from the metadata.")
    }
    
    dt_final
}

# ==========================================================================
# Loading Taxonomy
# ==========================================================================
 
#' Load and format a taxonomy reference table
#'
#' Reads a taxonomy TSV file, standardizes column names to lowercase, and extracts
#' the join column and taxonomic ranks based on the specified source type 
#' (`"reads"`, `"contigs"`, or `"kegg"`). Deduplicates entries on the key column
#' to ensure unique mapping.
#'
#' @param taxonomy_path Character. Path to the taxonomy TSV file.
#' @param source Character. Type of data source (`"reads"`, `"contigs"`, or `"kegg"`). 
#'   Determines the join key (`tax_id`, `contig_id`, or `ko`) and the expected 
#'   taxonomic ranks.
#'
#' @return A named \code{list} with three elements:
#'   \item{dt}{A \code{data.table} containing the unique taxonomy records.}
#'   \item{join_col}{Character. Name of the column used for merging/joining.}
#'   \item{tax_ranks}{Character vector. Names of the available taxonomic ranks present in the file.}
load_taxonomy_table <- function(taxonomy_path, source) {
  
  if (!file.exists(taxonomy_path)) {
    stop(sprintf("❌ ERROR: Taxonomy file not found: %s", taxonomy_path))
  }

  dt_taxo <- data.table::fread(taxonomy_path, showProgress = FALSE)
  data.table::setnames(dt_taxo, trimws(tolower(names(dt_taxo))))

  if (is.null(source) || is.na(source) || trimws(source) == "") {
  stop("❌ ERROR: 'source' argument is missing or empty in load_taxonomy_table().")
  }

  src <- tolower(trimws(source))

  if (grepl("reads", src)) {
    join_col  <- "tax_id"
    tax_ranks <- c("tax_id", "scientific_name", "domain", "kingdom",
                   "phylum", "class", "order", "family", "genus", "species")
  } else if (grepl("contigs", src)) {
    join_col  <- "contig_id"
    tax_ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
  } else if (grepl("kegg", src)) {
    join_col  <- "ko"
    tax_ranks <- c("ec_number", "level_1", "level_2", "level_3", "gene_description")
  } else {
    stop(sprintf("❌ ERROR: Unknown source [%s]. Expected 'reads', 'contigs', or 'kegg'.", source))
  }

  if (!join_col %in% names(dt_taxo)) {
    stop(sprintf(
      "❌ ERROR: Expected join column '%s' not found in taxonomy file (%s). Available columns: %s",
      join_col, taxonomy_path, paste(names(dt_taxo), collapse = ", ")
    ))
  }

  # Conversion sécurisée sans get()
  data.table::set(dt_taxo, j = join_col, value = as.character(dt_taxo[[join_col]]))

  if (!grepl("kegg", src) && anyDuplicated(dt_taxo[[join_col]])) {
    n_dup <- sum(duplicated(dt_taxo[[join_col]]))
    warning(sprintf(
      "⚠️ WARNING: %d duplicated '%s' found in taxonomy file — keeping first occurrence only.",
      n_dup, join_col
    ))
    dt_taxo <- unique(dt_taxo, by = join_col)
  }

  tax_ranks_present <- base::intersect(tax_ranks, names(dt_taxo))
  cols_to_keep      <- unique(c(join_col, tax_ranks_present))

  list(
    dt        = dt_taxo[, cols_to_keep, with = FALSE],
    join_col  = join_col,
    tax_ranks = tax_ranks_present
  )
}

annotate_with_taxonomy <- function(dt_results, taxo_ref, id_col_results = "Feature_ID") {
  
  if (is.null(dt_results) || !is.data.frame(dt_results) || NCOL(dt_results) == 0 || NROW(dt_results) == 0) {
    return(dt_results)
  } # enlever 

  # Auto-détection de la colonne ID dans le tableau de résultats
  if (is.null(id_col_results)) {
    possible_ids <- c("Feature_ID", "feature_id", taxo_ref$join_col, "kegg_id", "read_id", "contig_id", "rn")
    id_col_results <- intersect(possible_ids, names(dt_results))[1]
    
    if (is.na(id_col_results)) {
      stop(sprintf(
        "❌ ERROR: Impossible de trouver la colonne ID dans les résultats. Colonnes disponibles : %s",
        paste(names(dt_results), collapse = ", ")
      ))
    }
  } else if (!id_col_results %in% names(dt_results)) {
    stop(sprintf(
      "❌ ERROR: Column '%s' not found in results table. Available columns: %s",
      id_col_results, paste(names(dt_results), collapse = ", ")
    ))
  }
  # Conversion sécurisée sans get()
  data.table::set(dt_results, j = id_col_results, value = as.character(dt_results[[id_col_results]]))

  dt_annotated <- merge(
    dt_results, taxo_ref$dt,
    by.x = id_col_results, by.y = taxo_ref$join_col, all.x = TRUE
  )

  # Remplacement propre des NA par "Unclassified"
  for (col in taxo_ref$tax_ranks) {
    if (col %in% names(dt_annotated)) {
      na_idx <- which(is.na(dt_annotated[[col]]) | dt_annotated[[col]] == "")
      if (length(na_idx) > 0) {
        data.table::set(dt_annotated, i = na_idx, j = col, value = "Unclassified")
      }
    }
  }
  
  return(dt_annotated)
}

# ==========================================================================
# Loading Counts for DESeq2
# ==========================================================================

#' Build a DESeq2-ready count matrix from a list of per-sample TSV files.
#'
#' Matches each file to a sample_id via whole-token matching against
#' valid_ids (escaped regex, bounded by non-alphanumeric separators).
#' Unlike .match_sample_id(), this raises stop() (not warning()) on any
#' ambiguous match, since a wrong sample assignment here would silently
#' corrupt downstream DESeq2 results.
#'
#' @param data_paths Vector of paths to the TSV files.
#' @param valid_ids Character vector of valid sample_id values (typically
#'   meta_dt$sample_id).
#' @param id_cols Candidate feature-ID column names, tried in order.
#' @param count_col Name of the column holding raw counts.
#' @return A numeric count matrix (features x samples), NA replaced by 0.
build_deseq_count_matrix <- function(data_paths,
                                      valid_ids,
                                      id_cols = c("read_id", "contig_id", "kegg_id"),
                                      count_col = "read_mapped") {

  raw_list <- lapply(data_paths, function(f) {
    file_name <- basename(f)

    matched_sample <- valid_ids[vapply(valid_ids, function(sid) {
      sid_escaped <- .regex_escape(sid)
      pattern <- paste0("(^|[_.-])", sid_escaped, "([_.-]|$)")
      grepl(pattern, file_name, ignore.case = TRUE)
    }, FUN.VALUE = logical(1))]

    if (length(matched_sample) == 0 || is.na(matched_sample[1]) || matched_sample[1] == "") {
      stop(paste0(
        "\n❌ ERROR: No matching sample_id found in metadata for file: '", file_name, "'\n",
        "Please check if this sample is declared in your metadata.xlsx or verify the filename."
      ))
    }

    # Safeguard 1: Ambiguous matching detection
    if (length(matched_sample) > 1) {
      stop(paste0(
        "\n❌ ERROR: Multiple sample_ids matched for file '", file_name, "': ",
        paste(matched_sample, collapse = ", "), ". Please refine sample naming."
      ))
    }

    dt <- data.table::fread(f, showProgress = FALSE)
    if (nrow(dt) == 0) return(NULL)

    data.table::setnames(dt, tolower(names(dt)))

    current_id <- intersect(id_cols, names(dt))[1]
    if (is.na(current_id)) return(NULL)

    # Safeguard 2: Check for duplicated feature IDs to prevent cartesian explosion during merge
    if (anyDuplicated(dt[[current_id]])) {
      stop(paste0(
        "❌ ERROR: Duplicated feature IDs found in column '", current_id,
        "' for file: '", file_name, "'. Merging would duplicate rows and distort counts."
      ))
    }

    # Safeguard 3: Verify existence of the count column
    if (!count_col %in% names(dt)) {
      stop(paste0(
        "❌ ERROR: Column '", count_col, "' not found in file: '", file_name, "'. ",
        "Available columns: ", paste(names(dt), collapse = ", ")
      ))
    }

    dt[, .(
      feature_id = as.character(get(current_id)),
      sample_id  = matched_sample[1],
      count      = as.numeric(get(count_col))
    )]
  })

  raw_list <- Filter(Negate(is.null), raw_list)
  if (length(raw_list) == 0) {
    stop("🚨 Step error: raw_list is empty. No valid sample data tables were loaded for DESeq2 analysis.")
  }

  long_dt <- data.table::rbindlist(raw_list)
  rm(raw_list)

  count_data <- data.table::dcast(
    long_dt,
    feature_id ~ sample_id,
    value.var = "count",
    fill = 0
  )
  rm(long_dt)

  count_matrix <- as.matrix(count_data[, !"feature_id", with = FALSE])
  rownames(count_matrix) <- count_data$feature_id
  count_matrix <- round(count_matrix)

  n_na <- sum(is.na(count_matrix))
  if (n_na > 0) {
    warning("⚠️ ", n_na, " NA values detected after pivoting (likely due to a failure to parse numeric values) — replaced with 0.")
  }
  count_matrix[is.na(count_matrix)] <- 0

  count_matrix
}