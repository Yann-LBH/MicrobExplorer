################################################################################
# Project : "MicrobExplorer"
# Script: " utils loading metadata and data"
# Author: "Yann Le Bihan"
# Date: "2025-12-01"
# Link : https://github.com/Yann-LBH/MicrobExplorer
################################################################################

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
# Loading Read Counts for DESeq2
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
      grepl(paste0("(^|[^A-Za-z0-9_])", sid_escaped, "([^A-Za-z0-9_]|$)"), file_name)
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