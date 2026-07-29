test_that("get_GO_info() returns a term/description/ontology lookup table", {
  info <- get_GO_info()
  expect_s3_class(info, "tbl_df")
  expect_setequal(colnames(info), c("term", "description", "ontology"))
  expect_true("GO:0008150" %in% info$term)
  expect_true(all(c("BP", "CC", "MF") %in% unique(info$ontology)))
})

test_that("get_GO_levels() returns levels 1-3 for the three root terms' descendants", {
  levels <- get_GO_levels()
  expect_s3_class(levels, "tbl_df")
  expect_setequal(colnames(levels), c("GO_ID", "GO_lvl"))
  expect_setequal(
    levels$GO_lvl[
      levels$GO_ID %in% c("GO:0008150", "GO:0003674", "GO:0005575")
    ],
    1
  )
  expect_true(all(levels$GO_lvl %in% 1:3))
})
