# Load packages
if (!"org.Hs.eg.db" %in% installed.packages()) {
  BiocManager::install("org.Hs.eg.db")
}
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(enrichhelper)

# Get all human Entrez Gene IDs from the current reference annotation (org.Hs.eg.db)
human_entrez_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Make a mock DE result table with 20,000 genes.
# To produce a non-trivial (i.e., non-empty) enrichment result, inject real signal:
# genes annotated to a real GO:BP term ("transcription by RNA polymerase I", GO:0006360)
# are given low (significant) padj values and a positive lfc.
# All other genes get random padj/lfc values, i.e. no signal.
set.seed(4187)
signal_term <- "GO:0006360"
signal_genes <- AnnotationDbi::as.list(org.Hs.egGO2ALLEGS)[[signal_term]] |>
  unique()

DE_res <- tibble(
  contrast = "example_contrast",
  gene = human_entrez_ids[1:20000],
  padj = runif(20000),
  lfc = runif(20000, min = -20, max = 20)
) |>
  mutate(
    padj = ifelse(gene %in% signal_genes, runif(n(), 0, 0.01), padj),
    lfc = ifelse(gene %in% signal_genes, runif(n(), 5, 20), lfc)
  )

# Run GO enrichment analysis with the mock DE result table and the human GO mapping.
# NOTE: GO_ontology is ignored here since return_df = TRUE always runs all 3 ontologies
# (BP, MF, CC); GO_ontology only applies when return_df = FALSE.
GO_res <- run_ora(
  df = DE_res,
  contrast = "example_contrast",
  DE_direction = "up",
  ontology_type = "GO",
  OrgDb = "org.Hs.eg.db",
  simplify_terms = TRUE,
  return_df = TRUE
)

# Plot
GO_res |>
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
