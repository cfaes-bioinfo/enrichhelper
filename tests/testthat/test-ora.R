# Shared small term_map: 2 terms, overlapping gene membership
make_term_map <- function() {
  tibble::tibble(
    term = c(rep("TERM_A", 15), rep("TERM_B", 15)),
    gene = c(paste0("g", 1:15), paste0("g", 11:25)),
    description = c(
      rep("Term A description", 15),
      rep("Term B description", 15)
    )
  )
}

test_that("mean_lfc/median_lfc reflect only the requested contrast", {
  term_map <- make_term_map()
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]

  de_multi <- dplyr::bind_rows(
    tibble::tibble(
      contrast = "c1",
      gene = all_genes,
      lfc = ifelse(gene %in% term_a_genes, 10, 0.1),
      padj = ifelse(gene %in% term_a_genes, 0.001, 0.9)
    ),
    tibble::tibble(contrast = "c2", gene = all_genes, lfc = -10, padj = 0.9)
  )

  res <- run_ora(
    df = de_multi,
    term_map = term_map,
    contrast = "c1",
    DE_direction = "either",
    return_df = TRUE,
    verbose = FALSE
  )
  term_a_row <- res |> dplyr::filter(term == "TERM_A")

  expect_equal(nrow(term_a_row), 1)
  expect_equal(term_a_row$mean_lfc, 10)
})

test_that("a 'log2FoldChange'-named LFC column works like an 'lfc'-named one", {
  term_map <- make_term_map()
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]
  de_c1 <- tibble::tibble(
    contrast = "c1",
    gene = all_genes,
    lfc = ifelse(gene %in% term_a_genes, 10, 0.1),
    padj = ifelse(gene %in% term_a_genes, 0.001, 0.9)
  )
  de_dxs <- de_c1 |> dplyr::rename(log2FoldChange = lfc)

  res_lfc <- run_ora(
    df = de_c1,
    term_map = term_map,
    contrast = "c1",
    DE_direction = "either",
    return_df = TRUE,
    verbose = FALSE
  )
  res_dxs <- expect_no_error(
    run_ora(
      df = de_dxs,
      term_map = term_map,
      contrast = "c1",
      DE_direction = "either",
      return_df = TRUE,
      verbose = FALSE
    )
  )

  expect_equal(
    res_dxs |> dplyr::filter(term == "TERM_A") |> dplyr::pull(mean_lfc),
    res_lfc |> dplyr::filter(term == "TERM_A") |> dplyr::pull(mean_lfc)
  )
})

test_that("'redundant' is never NaN/error, even if qvalue::qvalue() fails", {
  # A term_map this small (2 terms tested) routinely makes qvalue::qvalue() fail --
  # that's exactly the scenario we want to exercise.
  term_map <- make_term_map()
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]
  de_c1 <- tibble::tibble(
    contrast = "c1",
    gene = all_genes,
    lfc = ifelse(gene %in% term_a_genes, 10, 0.1),
    padj = ifelse(gene %in% term_a_genes, 0.001, 0.9)
  )

  res <- run_ora(
    df = de_c1,
    term_map = term_map,
    contrast = "c1",
    DE_direction = "either",
    return_df = TRUE,
    verbose = FALSE
  )
  expect_equal(is.na(res$redundant), res$padj >= 0.05)
})

test_that("filter_no_descrip drops the whole term, not just its description", {
  term_map <- make_term_map()
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]
  de_c1 <- tibble::tibble(
    contrast = "c1",
    gene = all_genes,
    lfc = ifelse(gene %in% term_a_genes, 10, 0.1),
    padj = ifelse(gene %in% term_a_genes, 0.001, 0.9)
  )
  term_map_nodesc <- term_map |>
    dplyr::mutate(
      description = ifelse(term == "TERM_B", NA_character_, description)
    )

  res <- run_ora(
    df = de_c1,
    term_map = term_map_nodesc,
    contrast = "c1",
    DE_direction = "either",
    filter_no_descrip = TRUE,
    return_df = TRUE,
    verbose = FALSE
  )
  expect_false("TERM_B" %in% res$term)
})

test_that("return_df = FALSE still carries a 'redundant' column on every tested term", {
  term_map <- make_term_map()
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]
  de_c1 <- tibble::tibble(
    contrast = "c1",
    gene = all_genes,
    lfc = ifelse(gene %in% term_a_genes, 10, 0.1),
    padj = ifelse(gene %in% term_a_genes, 0.001, 0.9)
  )

  res <- run_ora(
    df = de_c1,
    term_map = term_map,
    contrast = "c1",
    DE_direction = "either",
    return_df = FALSE,
    verbose = FALSE
  )
  expect_true("redundant" %in% colnames(res@result))
  expect_equal(nrow(res@result), length(unique(term_map$term)))
})

test_that("focal_genes without a 'contrast' does not error, and DE_direction is NA", {
  term_map <- make_term_map()
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]

  res <- expect_no_error(
    run_ora(
      focal_genes = term_a_genes[1:5],
      term_map = term_map,
      return_df = TRUE,
      verbose = FALSE
    )
  )
  expect_true("contrast" %in% colnames(res))
  expect_true(all(is.na(res$contrast)))
  expect_true(all(is.na(res$DE_direction)))
})

test_that("mean_lfc ignores NA lfc values rather than propagating NA", {
  term_map <- make_term_map()
  all_genes <- unique(term_map$gene)
  term_a_genes <- term_map$gene[term_map$term == "TERM_A"]
  de_na <- tibble::tibble(
    contrast = "c1",
    gene = all_genes,
    lfc = ifelse(gene %in% term_a_genes, 10, 0.1),
    padj = ifelse(gene %in% term_a_genes, 0.001, 0.9)
  )
  de_na$lfc[de_na$gene == "g1"] <- NA_real_

  res <- run_ora(
    df = de_na,
    term_map = term_map,
    contrast = "c1",
    DE_direction = "either",
    return_df = TRUE,
    verbose = FALSE
  )
  expect_false(is.na(
    res |> dplyr::filter(term == "TERM_A") |> dplyr::pull(mean_lfc)
  ))
})

test_that("run_ora() requires exactly one of 'term_map'/'OrgDb'", {
  expect_snapshot(error = TRUE, run_ora(focal_genes = c("g1", "g2")))
  expect_snapshot(
    error = TRUE,
    run_ora(
      focal_genes = c("g1", "g2"),
      term_map = make_term_map(),
      OrgDb = "org.Hs.eg.db"
    )
  )
})

test_that("run_ora() requires 'organism' for KEGG and 'OrgDb' for GO", {
  expect_snapshot(
    error = TRUE,
    run_ora(
      focal_genes = c("g1", "g2"),
      ontology_type = "KEGG",
      OrgDb = "org.Hs.eg.db"
    )
  )
  expect_snapshot(
    error = TRUE,
    run_ora(
      focal_genes = c("g1", "g2"),
      ontology_type = "GO",
      organism = "hsa"
    )
  )
})
