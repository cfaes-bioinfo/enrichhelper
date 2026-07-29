# Load packages
library(tidyverse)
library(clusterProfiler)
library(enrichhelper)

# Define input files
DE_file <- "example/DE_stinkbug.tsv"
GO_map_file <- "example/GO-map_stinkbug.tsv"

# Read input files
DE_res <- read_tsv(DE_file, show_col_types = FALSE)
GO_map <- read_tsv(GO_map_file, show_col_types = FALSE)

# GO enrichment analysis
GO_res <- run_ora(
  df = DE_res,
  term_map = GO_map,
  contrast = "long_short",
  DE_direction = "up",
  simplify_terms = TRUE,
  return_df = TRUE
)

# Subset GO results since there are too many significant terms to plot all of them
GO_sig <- GO_res |>
  # Keep only significant, non-redundant terms (redundant: NA = not significant)
  filter(!is.na(redundant), !redundant)
GO_sel <- GO_sig |>
  # Only take the 10 most significant terms per ontology
  slice_min(padj, n = 10, with_ties = FALSE, by = "ontology") |>
  # Sorting here will be respected in the plot - sort by fold enrichment, the axis variable:
  arrange(fold_enrich)

# Cleveland dotplot
GO_sig |>
  cdotplot(
    x_var = "fold_enrich",
    fill_var = "padj_log",
    label_var = "n_focal_in_cat",
    facet_var1 = "ontology",
    facet_to_columns = FALSE,
    add_term_id = FALSE,
    label_chars = 60,
    ylab_size = 10,
    scico_palette = "batlow"
  )
ggsave("example/cdotplot_stinkbug.png", width = 8, height = 6, dpi = 300)

# Cleveland dotplot of just the top-10-per-ontology subset
GO_sel |>
  cdotplot(
    x_var = "fold_enrich",
    fill_var = "padj_log",
    label_var = "n_focal_in_cat",
    facet_var1 = "ontology",
    facet_to_columns = FALSE,
    add_term_id = FALSE,
    label_chars = 60,
    ylab_size = 10,
    scico_palette = "batlow"
  )
