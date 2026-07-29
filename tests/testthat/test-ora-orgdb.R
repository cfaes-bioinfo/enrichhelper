test_that("return_df = TRUE runs all 3 GO ontologies regardless of GO_ontology", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  suppressPackageStartupMessages(library(org.Hs.eg.db))

  set.seed(2)
  entrez_ids <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")[1:3000]
  # A targeted select() query for just the one signal term, rather than
  # AnnotationDbi::as.list(org.Hs.egGO2ALLEGS) -- the latter materializes the full
  # gene-to-all-GO-ancestors mapping for every human gene (~1.5GB+) just to look up a
  # single term, which is unnecessarily memory-heavy for a unit test.
  signal_genes <- unique(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keytype = "GOALL",
      keys = "GO:0006360",
      columns = "ENTREZID"
    )$ENTREZID
  )
  de_human <- tibble::tibble(
    contrast = "c1",
    gene = entrez_ids,
    padj = stats::runif(length(entrez_ids)),
    lfc = stats::runif(length(entrez_ids), -5, 5)
  ) |>
    dplyr::mutate(
      padj = ifelse(
        gene %in% signal_genes,
        stats::runif(dplyr::n(), 0, 0.01),
        padj
      ),
      lfc = ifelse(gene %in% signal_genes, stats::runif(dplyr::n(), 5, 10), lfc)
    )

  res <- run_ora(
    df = de_human,
    contrast = "c1",
    DE_direction = "up",
    ontology_type = "GO",
    OrgDb = "org.Hs.eg.db",
    GO_ontology = "BP", # should be ignored -- all 3 ontologies expected regardless
    return_df = TRUE,
    verbose = FALSE
  )
  expect_setequal(unique(res$ontology), c("BP", "MF", "CC"))
})
