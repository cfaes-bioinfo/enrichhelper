test_that("run_gsea() finds an injected-signal term and returns a tidy df", {
  set.seed(42)
  term_map <- tibble::tibble(
    term = c(rep("TERM_A", 30), rep("TERM_B", 30)),
    gene = c(paste0("g", 1:30), paste0("g", 51:80)),
    description = c(
      rep("Term A description", 30),
      rep("Term B description", 30)
    )
  )
  other_genes <- paste0("g", setdiff(1:150, c(1:30, 51:80)))
  de_df <- tibble::tibble(
    contrast = "c1",
    gene = c(paste0("g", 1:30), other_genes, paste0("g", 51:80)),
    lfc = c(
      stats::runif(30, 5, 10),
      stats::runif(length(other_genes), -2, 2),
      stats::runif(30, -10, -5)
    ),
    padj = stats::runif(150)
  )

  # fgsea emits benign p-value-precision warnings at this small a scale (permutation
  # counts are low with only 2 pathways/150 genes) -- suppress them here.
  res <- suppressWarnings(run_gsea(
    df = de_df,
    contrast = "c1",
    term_map = term_map,
    return_df = TRUE
  ))

  expect_s3_class(res, "tbl_df")
  expect_true(all(
    c("term", "padj", "sig", "description", "mean_lfc", "median_lfc") %in%
      colnames(res)
  ))
  term_a <- res |> dplyr::filter(term == "TERM_A")
  expect_equal(nrow(term_a), 1)
  expect_true(term_a$sig)
  expect_gt(term_a$mean_lfc, 0)
})
