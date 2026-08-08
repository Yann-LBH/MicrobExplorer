################################################################################
# Project : "MicrobExplorer"
# Script  : "Utils resolve title for plots"
# Author  : "Yann Le Bihan"
# Date    : "2025-12-01"
# Link    : https://github.com/Yann-LBH/MicrobExplorer
#
################################################################################
suppressPackageStartupMessages({
  library(jsonlite)
  library(glue)
})

# Il va directement piocher dans snakemake
JSON  <- as.character(snakemake@input[["title_resolver"]])[1]
.lang <- as.character(snakemake@params[["language"]])[1]

if (!file.exists(JSON)) {
  stop(sprintf("Fichier de traduction JSON introuvable : %s", JSON))
}

.title_dict <- jsonlite::fromJSON(JSON, simplifyDataFrame = TRUE)

t <- function(key) {
  row <- .title_dict $translation[.title_dict $translation$key == key, ]
  if (nrow(row) == 0) stop(sprintf("Clé de traduction '%s' introuvable.", key))
  val <- row[[.lang]]
  if (is.null(val) || is.na(val)) stop(sprintf("Langue '%s' absente pour '%s'.", .lang, key))
  return(val)
}

resolve_text <- function(key, ...) {
  template <- t(key)
  glue::glue(template, ..., .envir = parent.frame())
}