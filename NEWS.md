# enrichhelper 0.0.0.9004

* `run_gsea()` gained support for KEGG GSEA of model organisms (`ontology_type
  = "KEGG"`) via `clusterProfiler::gseKEGG()`, using a new `organism` argument
  (KEGG species code, e.g. `"hsa"`) instead of a `term_map`, mirroring the
  `run_ora()` KEGG support added in 0.0.0.9003. `simplify_terms` is not
  supported for KEGG, and an unsupported `keyType` (e.g. `"ENSEMBL"`) now
  errors immediately with a clear message. `run_gsea()` also gained the same
  upfront `term_map`/`OrgDb`/`organism`/`ontology_type` validation as
  `run_ora()`, rather than failing later with a less clear error.

# enrichhelper 0.0.0.9003

* `run_ora()` gained support for KEGG enrichment of model organisms
  (`ontology_type = "KEGG"`) via `clusterProfiler::enrichKEGG()`, using a new
  `organism` argument (KEGG species code, e.g. `"hsa"`) instead of a
  `term_map`. `simplify_terms` is not supported for KEGG (no GO-like term
  structure). A `keyType` unsupported by KEGG's own API (e.g. `"ENSEMBL"`)
  now errors immediately with a clear message, rather than surfacing a
  cryptic error from deep inside `enrichKEGG()`.

# enrichhelper 0.0.0.9002

* `run_gsea()` gained a `simplify_terms` (+ `simplify_cutoff`) argument, mirroring
  `run_ora()`: uses `clusterProfiler::simplify()` to flag redundant/similar GO terms
  among the significant ones via a `redundant` column.
* `run_gsea()` gained an `n_perm` argument (`clusterProfiler`'s `nPermSimple`) to
  control the number of GSEA permutations, with a higher default (10,000) than
  `clusterProfiler`'s own default.

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
