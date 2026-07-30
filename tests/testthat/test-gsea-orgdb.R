test_that("run_gsea() runs KEGG GSEA for a model organism (ontology_type = 'KEGG')", {
  testthat::skip_if_not_installed("org.Hs.eg.db")
  testthat::skip_if_offline()
  suppressPackageStartupMessages(library(org.Hs.eg.db))

  set.seed(3)
  entrez_ids <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID")[1:3000]
  # Signal genes: members of the Glycolysis / Gluconeogenesis KEGG pathway (hsa00010),
  # via org.Hs.eg.db's own (older, but sufficiently overlapping) PATH mapping -- gseKEGG()
  # itself downloads the current gene set from KEGG at call time.
  signal_genes <- unique(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keytype = "PATH",
      keys = "00010",
      columns = "ENTREZID"
    )$ENTREZID
  )
  de_human <- tibble::tibble(
    contrast = "c1",
    gene = entrez_ids,
    padj = stats::runif(length(entrez_ids)),
    lfc = ifelse(
      entrez_ids %in% signal_genes,
      stats::runif(length(entrez_ids), 5, 10),
      stats::runif(length(entrez_ids), -2, 2)
    )
  )

  res <- suppressWarnings(run_gsea(
    df = de_human,
    contrast = "c1",
    ontology_type = "KEGG",
    organism = "hsa",
    return_df = TRUE
  ))
  top_term <- res |> dplyr::filter(padj == min(padj)) |> dplyr::pull(term)
  expect_identical(top_term, "hsa00010")
})
