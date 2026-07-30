# Shared fixture: 2 terms with an injected LFC signal, used across tests below
make_gsea_fixture <- function() {
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
  list(term_map = term_map, de_df = de_df)
}

test_that("run_gsea() finds an injected-signal term and returns a tidy df", {
  set.seed(42)
  fixture <- make_gsea_fixture()

  # fgsea emits benign p-value-precision warnings at this small a scale (permutation
  # counts are low with only 2 pathways/150 genes) -- suppress them here.
  res <- suppressWarnings(run_gsea(
    df = fixture$de_df,
    contrast = "c1",
    term_map = fixture$term_map,
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

test_that("run_gsea() folds significance into 'redundant' when simplify_terms = FALSE", {
  set.seed(42)
  fixture <- make_gsea_fixture()

  res <- suppressWarnings(run_gsea(
    df = fixture$de_df,
    contrast = "c1",
    term_map = fixture$term_map,
    return_df = TRUE
  ))

  expect_true("redundant" %in% colnames(res))
  # NA = not significant, FALSE = significant (never TRUE without simplify_terms = TRUE)
  expect_equal(is.na(res$redundant), res$padj >= 0.05)
  expect_false(any(res$redundant, na.rm = TRUE))
})

test_that("run_gsea() with return_df = FALSE drops non-significant terms and carries 'redundant'", {
  set.seed(42)
  fixture <- make_gsea_fixture()

  res <- suppressWarnings(run_gsea(
    df = fixture$de_df,
    contrast = "c1",
    term_map = fixture$term_map,
    return_df = FALSE
  ))

  expect_true("redundant" %in% colnames(res@result))
  expect_true(all(res@result$p.adjust < 0.05))
  expect_false(any(res@result$redundant))
})

test_that("run_gsea() warns and skips simplify when term_map has no 'ontology' column", {
  set.seed(42)
  fixture <- make_gsea_fixture()

  # fgsea can additionally emit benign p-value-precision warnings at this small a scale
  # -- muffle those so only the 'ontology' warning under test reaches expect_warning().
  expect_warning(
    withCallingHandlers(
      res <- run_gsea(
        df = fixture$de_df,
        contrast = "c1",
        term_map = fixture$term_map,
        simplify_terms = TRUE,
        return_df = TRUE
      ),
      warning = function(w) {
        if (!grepl("ontology", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ),
    "ontology"
  )
  expect_true("redundant" %in% colnames(res))
  expect_false(any(res$redundant, na.rm = TRUE))
})

test_that("run_gsea() requires exactly one of 'term_map'/'OrgDb'/'organism'", {
  fixture <- make_gsea_fixture()
  expect_snapshot(
    error = TRUE,
    run_gsea(df = fixture$de_df, contrast = "c1")
  )
  expect_snapshot(
    error = TRUE,
    run_gsea(
      df = fixture$de_df,
      contrast = "c1",
      term_map = fixture$term_map,
      OrgDb = "org.Hs.eg.db"
    )
  )
})

test_that("run_gsea() requires 'organism' for KEGG and 'OrgDb' for GO", {
  fixture <- make_gsea_fixture()
  expect_snapshot(
    error = TRUE,
    run_gsea(
      df = fixture$de_df,
      contrast = "c1",
      ontology_type = "KEGG",
      OrgDb = "org.Hs.eg.db"
    )
  )
  expect_snapshot(
    error = TRUE,
    run_gsea(
      df = fixture$de_df,
      contrast = "c1",
      ontology_type = "GO",
      organism = "hsa"
    )
  )
})

test_that("run_gsea() rejects a KEGG keyType unsupported by KEGG's own API", {
  fixture <- make_gsea_fixture()
  expect_snapshot(
    error = TRUE,
    run_gsea(
      df = fixture$de_df,
      contrast = "c1",
      ontology_type = "KEGG",
      organism = "hsa",
      keyType = "ENSEMBL"
    )
  )
})
