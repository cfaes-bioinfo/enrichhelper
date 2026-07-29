# Load packages
library(tidyverse)
library(DESeq2)
library(enrichhelper)

# Define input files
dds_file <- "example/deseq_object.rds"
GO_map_file <- "example/GO-map_stinkbug.tsv"

# Read input files
dds <- readRDS(dds_file)
GO_map <- read_tsv(GO_map_file, show_col_types = FALSE)

# DE analysis for adults only, by daylength
dds <- dds[, dds$age == "adult"]
design(dds) <- formula(~daylength)
dds <- DESeq(dds)

# Extract results
res_deseq <- results(dds, contrast = c("daylength", "long", "short"))

# run_ora() accepts a raw DESeq2 results() object directly -- it's converted
# internally (gene IDs pulled from rownames, log2FoldChange/baseMean/pvalue
# renamed to lfc/mean/p) and tagged with the 'contrast' argument below.
GO_res <- run_ora(
  df = res_deseq,
  term_map = GO_map,
  contrast = "long_short",
  DE_direction = "up",
  simplify_terms = TRUE,
  return_df = TRUE
)
