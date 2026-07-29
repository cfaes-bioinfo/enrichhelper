make_enrich_df <- function() {
  tibble::tibble(
    contrast = "c1",
    DE_direction = "up",
    term = c("GO:0000001", "GO:0000002", "GO:0000003"),
    description = c(
      "term one description",
      "term two description",
      NA_character_
    ),
    padj = c(0.001, 0.01, 0.02),
    n_focal_in_cat = c(5, 3, 2),
    median_lfc = c(1.5, -2, 0.5),
    mean_lfc = c(1.4, -1.9, 0.6),
    fold_enrich = c(2, 3, 1.5)
  )
}

test_that("cdotplot() returns a ggplot object for a minimal valid df", {
  p <- cdotplot(make_enrich_df())
  expect_s3_class(p, "ggplot")
})

test_that("cdotplot() errors with a clear message when a required column is missing", {
  df <- make_enrich_df() |> dplyr::select(-term)
  expect_snapshot(error = TRUE, cdotplot(df))
})

test_that("cdotplot() errors when DE_dirs is given but DE_direction is missing", {
  df <- make_enrich_df() |> dplyr::select(-DE_direction)
  expect_snapshot(error = TRUE, cdotplot(df, DE_dirs = "up"))
})

test_that("cdotplot() errors when x_var/fill_var/facet_var names a missing column", {
  expect_snapshot(
    error = TRUE,
    cdotplot(make_enrich_df(), facet_var1 = "ontology")
  )
})

test_that("cdotplot() errors when nothing is left to plot after filtering", {
  df <- make_enrich_df() |> dplyr::filter(FALSE)
  expect_snapshot(error = TRUE, cdotplot(df))
})
