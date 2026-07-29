# A minimal stand-in for a real DESeq2 results() object: prep_deseq_df() (run_ora()'s
# DESeqResults handler) only needs the object to (a) satisfy inherits(df, "DESeqResults")
# and (b) convert via as.data.frame() to a data.frame with baseMean/log2FoldChange/
# lfcSE/stat/pvalue/padj columns and gene IDs as rownames -- exactly what a real
# DESeqResults gives via its as.data.frame() method. Building this by hand (instead of
# fitting a real DESeqDataSet) keeps this test fast and avoids pulling in DESeq2's large
# dependency tree (S4Vectors/IRanges/GenomicRanges/BiocParallel/...) just for a plumbing
# test that isn't exercising DESeq2's own statistics.
make_fake_deseq_results <- function(
  gene,
  baseMean,
  log2FoldChange,
  pvalue,
  padj
) {
  df <- data.frame(
    baseMean = baseMean,
    log2FoldChange = log2FoldChange,
    lfcSE = 0.2,
    stat = log2FoldChange / 0.2,
    pvalue = pvalue,
    padj = padj,
    row.names = gene
  )
  class(df) <- c("DESeqResults", "data.frame")
  df
}

test_that("a DESeqResults-like object is accepted and tagged with 'contrast'", {
  term_map <- tibble::tibble(
    term = c(rep("TERM_A", 15), rep("TERM_B", 15)),
    gene = c(paste0("g", 1:15), paste0("g", 11:25)),
    description = c(
      rep("Term A description", 15),
      rep("Term B description", 15)
    )
  )
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]
  is_term_a <- all_genes %in% term_a_genes

  res_deseq <- make_fake_deseq_results(
    gene = all_genes,
    baseMean = 50,
    log2FoldChange = ifelse(is_term_a, 3, 0.1),
    pvalue = ifelse(is_term_a, 1e-10, 0.8),
    padj = ifelse(is_term_a, 1e-9, 0.9)
  )
  expect_true(inherits(res_deseq, "DESeqResults"))

  res <- expect_no_error(
    run_ora(
      df = res_deseq,
      term_map = term_map,
      contrast = "trt_vs_ctrl",
      DE_direction = "either",
      return_df = TRUE,
      verbose = FALSE
    )
  )
  expect_true(all(res$contrast == "trt_vs_ctrl"))
  expect_false(res |> dplyr::filter(term == "TERM_A") |> dplyr::pull(redundant))
})
