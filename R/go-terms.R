#' Get descriptions and ontologies for all GO terms
#'
#' Returns a lookup table with the description and ontology (`BP`, `MF`, or
#' `CC`) for every GO term in the installed `GO.db` annotation package.
#'
#' @return A tibble with columns `term`, `description`, and `ontology`.
#' @export
#' @examples
#' \donttest{
#' GO_info <- get_GO_info()
#' }
get_GO_info <- function() {
  AnnotationDbi::select(
    GO.db::GO.db,
    columns = c("GOID", "TERM", "ONTOLOGY"),
    keys = AnnotationDbi::keys(GO.db::GO.db, keytype = "GOID"),
    keytype = "GOID"
  ) |>
    dplyr::rename(
      term = GOID,
      description = TERM,
      ontology = ONTOLOGY
    ) |>
    tibble::as_tibble()
}

# Fetch the direct children of 'go_ids' across all three GO sub-ontologies
# (BP/CC/MF), dropping any IDs with no children. Used by get_GO_levels().
# NOTE: must call AnnotationDbi::mget() (an S4 generic dispatching on the
# GOBPCHILDREN/... Bimap objects), not base::mget() -- the latter doesn't know
# how to handle those objects, and (unlike in a plain script) attaching GO.db/
# AnnotationDbi to the search path doesn't make base::mget() dispatch correctly
# here, because package code doesn't fall back to the search path for unqualified
# lookups the way top-level script code does.
get_GO_children <- function(go_ids) {
  bp <- unname(unlist(AnnotationDbi::mget(
    go_ids,
    GO.db::GOBPCHILDREN,
    ifnotfound = NA
  )))
  cc <- unname(unlist(AnnotationDbi::mget(
    go_ids,
    GO.db::GOCCCHILDREN,
    ifnotfound = NA
  )))
  mf <- unname(unlist(AnnotationDbi::mget(
    go_ids,
    GO.db::GOMFCHILDREN,
    ifnotfound = NA
  )))

  all_children <- c(bp, cc, mf)
  unique(all_children[!is.na(all_children)])
}

#' Get GO levels for the most generic GO terms
#'
#' Builds a lookup table with the GO "level" (distance from the root term of
#' its sub-ontology) for the three root terms (level 1) and their level-2 and
#' level-3 descendants across all three GO sub-ontologies (`BP`, `MF`, `CC`).
#' More specific (deeper) terms are not currently included.
#'
#' @return A tibble with columns `GO_ID` and `GO_lvl`.
#' @export
#' @examples
#' \donttest{
#' GO_levels <- get_GO_levels()
#' }
get_GO_levels <- function() {
  # The three root terms
  lvl1_ids <- c("GO:0008150", "GO:0003674", "GO:0005575")
  lvl2_ids <- get_GO_children(lvl1_ids)
  lvl3_ids <- get_GO_children(lvl2_ids)

  tibble::tibble(GO_ID = unique(c(lvl1_ids, lvl2_ids, lvl3_ids))) |>
    dplyr::mutate(
      GO_lvl = dplyr::case_when(
        GO_ID %in% lvl1_ids ~ 1,
        GO_ID %in% lvl2_ids ~ 2,
        GO_ID %in% lvl3_ids ~ 3,
        TRUE ~ NA_real_
      )
    )
}
