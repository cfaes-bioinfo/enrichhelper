#' @keywords internal
#' @importFrom dplyr any_of arrange case_when desc distinct filter
#' @importFrom dplyr join_by left_join mutate pull relocate rename right_join
#' @importFrom dplyr select summarize
#' @importFrom tibble tibble as_tibble rownames_to_column
#' @importFrom tidyr separate_longer_delim separate_wider_delim
#' @importFrom purrr map_dfr
#' @importFrom stringr str_trunc
#' @importFrom ggplot2 element_blank element_text expansion facet_grid
#' @importFrom ggplot2 geom_hline labs margin scale_color_brewer
#' @importFrom ggplot2 scale_color_viridis_c scale_x_discrete
#' @importFrom ggplot2 scale_y_continuous theme theme_bw vars
#' @importFrom methods is
#' @importFrom rlang .data
"_PACKAGE"

# NSE column names used unquoted inside dplyr/tidyr/ggplot2 verbs throughout
# the package -- declared here so R CMD check doesn't flag them as undefined
# global variables (see https://r-pkgs.org/data.html#sec-data-globalvars).
utils::globalVariables(c(
  "GO_ID",
  "GOID",
  "TERM",
  "p.adjust",
  "qvalue",
  "Count",
  "ID",
  "GeneRatio",
  "BgRatio",
  "Description",
  "geneID",
  "gene_ids",
  "gene",
  "lfc",
  "term",
  "ontology",
  "GeneRatio_1",
  "GeneRatio_2",
  "BgRatio_1",
  "BgRatio_2",
  "n_focal_in_cat",
  "n_focal",
  "n_cat",
  "n_total",
  "padj",
  "isDE",
  "log2FoldChange",
  "baseMean",
  "pvalue",
  "ONTOLOGY",
  "core_enrichment",
  "sig",
  "description",
  "contrast",
  "DE_direction",
  "redundant"
))
