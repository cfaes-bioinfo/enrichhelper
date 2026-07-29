# enrichhelper 0.0.0.9001

* Fixed `remotes::install_github()` failing with a git error for `GO.db`
  (an annotation data package not hosted on Bioconductor's git server).
  Replaced the `Remotes: bioc::...` entries with a `biocViews` field, which
  lets `remotes`/`BiocManager` resolve Bioconductor dependencies from the
  standard package repositories instead.

# enrichhelper 0.0.0.9000

* Initial release, providing `run_ora()` and `run_gsea()` (wrappers around
  `clusterProfiler` for over-representation and gene set enrichment
  analysis), `cdotplot()` for visualizing enrichment results, `cp_to_df()`
  for converting a `clusterProfiler` result to a tidy tibble, and
  `get_GO_info()`/`get_GO_levels()` for GO term lookup tables.
