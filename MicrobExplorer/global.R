# Libraries
library(shiny)
library(bslib)
library(bsicons)
library(dplyr)
library(arrow)
#BiocManager::install("phyloseq")
#library(phyloseq)
library(ggplot2)
library(dplyr)
library(jsonlite)

# Source all module files automatically
module_files <- list.files("modules", full.names = TRUE, pattern = "\\.R$")
sapply(module_files, source)

# Load heavy input data once at startup
#parquet_files <- list.files(full.names = TRUE, pattern = "\\.parquet$")
# df_deseq <- read_parquet("data/deseq_results.parquet")
# ps_data  <- readRDS("data/phyloseq_object.rds")