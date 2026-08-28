#' Gene-length-aware over-representation analysis (a `goseq` reimplementation)
#'
#' A length-bias-corrected alternative to [run_ora()], for use when genes are
#' more likely to enter the focal set the longer they are (e.g. focal genes
#' called via coordinate overlap with a genomic feature). Reimplements the
#' method of the Bioconductor `goseq` package directly rather than depending
#' on it, because `goseq` needs `GenomicFeatures` -> `rtracklayer` -> `XML`,
#' and `XML` fails to install/load in some environments. The two real
#' dependencies of the method, `BiasedUrn` (the Wallenius distribution) and
#' `mgcv` (the spline), are used instead.
#'
#' `clusterProfiler::enricher()` (and so [run_ora()]) assumes every gene had
#' an equal chance of being selected into the focal set. When that is not
#' true -- for example when focal genes are called via coordinate overlap
#' with a genomic feature, so long genes are more likely to be selected --
#' the count of focal genes in a category is no longer central hypergeometric
#' under the null. `run_goseq()` corrects for this the way `goseq` does:
#' \enumerate{
#'   \item Fit a probability weighting function (PWF), P(focal) vs.
#'     `bias_data` (e.g. gene length), via a binomial GAM on
#'     `log(bias_data)` (`mgcv::gam()`). This differs from `goseq`'s own
#'     rolling-mean-plus-monotonic-cubic-spline PWF: the GAM keeps fitted
#'     values in (0, 1) by construction and does not impose monotonicity as
#'     an assumption. Set `pwf_k` to change the spline's flexibility.
#'   \item Per category, set odds `w` = mean(PWF inside) / mean(PWF outside).
#'   \item Get the p-value from the Wallenius noncentral hypergeometric
#'     (`BiasedUrn::pWNCHypergeo()`) with that odds.
#' }
#' A category full of long genes therefore gets a null that already expects a
#' high focal fraction, and only an excess beyond that counts as enrichment.
#'
#' Other differences from the `goseq` package: there are no q-values, so
#' `redundant` is `NA` on `padj` alone ([run_ora()] also requires `qvalue <
#' q_enrich`); and categories with no focal gene are dropped before BH
#' (`drop_empty_cats = TRUE`) to match `clusterProfiler`/[run_ora()] -- the
#' `goseq` package keeps them, which enlarges the BH denominator and makes it
#' more conservative (set `drop_empty_cats = FALSE` to match `goseq` instead).
#'
#' The returned tibble carries [run_ora()]'s columns, in [run_ora()]'s order,
#' so it is a drop-in for downstream code. goseq-specific columns
#' (`p_uncorrected`, `odds`, `method`) are appended at the end.
#'
#' @param df DE results df with at least columns `gene`, `contrast`,
#'   `lfc`/`log2FoldChange`, `padj`. Alternative to `focal_genes`; same rules
#'   as [run_ora()] for deriving the focal (DE) gene set.
#' @param contrast DE comparison as specified in the `contrast` column in `df`.
#' @param DE_direction DE direction: `"either"` (all DE genes), `"up"` (lfc >
#'   0), or `"down"` (lfc < 0). Ignored when `focal_genes` is supplied.
#' @param focal_genes Alternative to `df`: a pre-built vector of focal gene IDs.
#' @param term_map Term-to-gene map with columns: 1: term, 2: gene, and
#'   optionally 3: description, 4: ontology. Required.
#' @param universe Background universe of genes (vector of gene IDs). Always
#'   intersected with `term_map$gene` and with `names(bias_data)`; defaults
#'   to every gene in `term_map`.
#' @param bias_data Named numeric vector, gene ID -> bias covariate
#'   (typically gene length). Required. Universe genes missing a value are
#'   dropped (with a message) before the PWF is fit.
#' @param method `"Wallenius"` (length-corrected p-value) or
#'   `"Hypergeometric"` (plain hypergeometric p-value, i.e. `p_uncorrected`
#'   -- useful as a like-for-like uncorrected baseline).
#' @param p_enrich Adjusted p-value threshold for significance (used for the
#'   `redundant` column and the verbose enriched-term count).
#' @param p_DE Adjusted p-value threshold used to derive `isDE` from `df`
#'   when `df` has no `isDE` column.
#' @param min_DE_in_cat Min. number of focal genes in a category for
#'   significance (`redundant` column); also the minimum number of focal
#'   genes overall for the analysis to run at all.
#' @param min_cat_size Min. number of genes in a category.
#' @param max_cat_size Max. number of genes in a category.
#' @param drop_empty_cats Drop categories with zero focal genes before BH
#'   correction (matches `clusterProfiler`/[run_ora()]; set `FALSE` to match
#'   the `goseq` package, which keeps them and so is more conservative).
#' @param simplify_terms Whether to use `clusterProfiler`'s `simplify()`
#'   (via its internal `simplify_internal()`) to flag redundant/similar GO
#'   terms among the significant ones. Requires an `ontology` column in
#'   `term_map`.
#' @param simplify_cutoff Simplify similarity cutoff.
#' @param pwf_k Basis dimension (`k`) of the PWF's `mgcv::gam()` smooth term;
#'   higher allows a more flexible P(focal) vs. `log(bias_data)` curve.
#' @param plot_pwf Plot the fitted PWF (P(focal) vs. `bias_data`, log-x) to
#'   the current graphics device.
#' @param allow_dups Allow a gene ID to be present multiple times in
#'   `focal_genes`/the derived focal gene list.
#' @param verbose Print a one-line progress message (contrast/DE direction/
#'   gene counts, and the enriched-term count -- with the non-redundant count
#'   alongside it when `simplify_terms = TRUE`).
#'
#' @return A tibble with one row per tested term, or `NULL` when there were
#'   too few focal genes (`< min_DE_in_cat`) or no categories passed
#'   `min_cat_size`/`max_cat_size`.
#' @export
run_goseq <- function(
  df = NULL,
  contrast = NULL,
  DE_direction = "either",
  focal_genes = NULL,
  term_map = NULL,
  universe = NULL,
  bias_data = NULL,
  method = c("Wallenius", "Hypergeometric"),
  p_enrich = 0.05,
  p_DE = 0.05,
  min_DE_in_cat = 2,
  min_cat_size = 5,
  max_cat_size = 500,
  drop_empty_cats = TRUE,
  simplify_terms = FALSE,
  simplify_cutoff = 0.7,
  pwf_k = 6,
  plot_pwf = FALSE,
  allow_dups = FALSE,
  verbose = TRUE
) {
  method <- match.arg(method)
  stopifnot(
    "Supply a 'term_map'" = !is.null(term_map),
    "Supply 'bias_data' as a named numeric vector (gene -> gene length)" =
      !is.null(bias_data) && !is.null(names(bias_data))
  )
  fcontrast <- if (is.null(contrast)) NA_character_ else contrast

  # -- Focal genes ------------------------------------------------------------
  # Same rules as enrichhelper:::prep_ora_genes(): derive isDE from padj when
  # absent, subset to the contrast, then to the requested LFC direction.
  if (is.null(focal_genes)) {
    stopifnot("Supply either 'focal_genes' or 'df'" = !is.null(df))
    if ("log2FoldChange" %in% colnames(df) && !"lfc" %in% colnames(df)) {
      df <- dplyr::rename(df, lfc = log2FoldChange)
    }
    if (!"isDE" %in% colnames(df)) {
      df <- dplyr::mutate(df, isDE = !is.na(padj) & padj < p_DE)
    }
    df_c <- dplyr::filter(df, contrast == fcontrast)
    stopifnot("No rows left after contrast filtering" = nrow(df_c) > 0)
    if (DE_direction == "up") df_c <- dplyr::filter(df_c, lfc > 0)
    if (DE_direction == "down") df_c <- dplyr::filter(df_c, lfc < 0)
    focal_genes <- df_c |> dplyr::filter(isDE) |> dplyr::pull(gene)
  } else {
    DE_direction <- NA_character_
  }
  if (any(duplicated(focal_genes))) {
    if (allow_dups) {
      focal_genes <- unique(focal_genes)
    } else {
      stop("Duplicated gene IDs in focal_genes -- multiple contrasts in 'df'?")
    }
  }

  # -- Universe ---------------------------------------------------------------
  # Always intersected with the term map and with bias_data: a gene with no
  # length cannot be weighted, and including it would silently distort the PWF.
  univ <- if (is.null(universe)) unique(term_map$gene) else unique(universe)
  univ <- univ[univ %in% term_map$gene]
  n_dropped <- sum(!univ %in% names(bias_data))
  if (n_dropped > 0 && verbose) {
    message("Note: dropped ", n_dropped, " universe genes with no bias_data value")
  }
  univ <- univ[univ %in% names(bias_data)]
  focal <- unique(focal_genes[focal_genes %in% univ])

  if (length(focal) < min_DE_in_cat) {
    if (verbose) {
      cat("Contrast:", fcontrast, "// DE direction:", DE_direction,
          "// too few focal genes (", length(focal), ") -- returning NULL\n")
    }
    return(NULL)
  }

  # -- 1. Probability weighting function --------------------------------------
  pwf_df <- data.frame(
    gene = univ,
    len = as.numeric(bias_data[univ]),
    focal = as.integer(univ %in% focal)
  )
  gam_fit <- mgcv::gam(
    focal ~ s(log(len), k = pwf_k),
    family = stats::binomial(),
    data = pwf_df
  )
  pwf_df$pwf <- as.numeric(stats::fitted(gam_fit))
  # Guard against 0/1 fitted values, which would make the odds 0 or infinite
  pwf_df$pwf <- pmin(pmax(pwf_df$pwf, 1e-6), 1 - 1e-6)
  pwf <- stats::setNames(pwf_df$pwf, pwf_df$gene)

  if (plot_pwf) {
    o <- order(pwf_df$len)
    plot(pwf_df$len[o], pwf_df$pwf[o], type = "l", log = "x", ylim = c(0, 1),
         xlab = "Gene length (bp)", ylab = "P(focal)",
         main = paste("PWF:", fcontrast, DE_direction))
    rug(pwf_df$len[pwf_df$focal == 1], col = "red")
  }

  # -- Category membership ----------------------------------------------------
  tg <- term_map |>
    dplyr::filter(gene %in% univ) |>
    dplyr::distinct(term, gene)
  sets <- split(tg$gene, tg$term)
  sizes <- lengths(sets)
  sets <- sets[sizes >= min_cat_size & sizes <= max_cat_size]
  if (length(sets) == 0) return(NULL)

  N <- length(univ)
  n_focal <- length(focal)
  sum_pwf_all <- sum(pwf)

  # -- 2./3. Per-category odds and p-value ------------------------------------
  per_cat <- function(genes) {
    n_cat <- length(genes)
    hits <- genes[genes %in% focal]
    x <- length(hits)
    s_in <- sum(pwf[genes])
    # Wallenius odds: average selection weight inside vs. outside the category
    odds <- (s_in / n_cat) / ((sum_pwf_all - s_in) / (N - n_cat))
    p_hyper <- stats::phyper(x - 1, n_cat, N - n_cat, n_focal, lower.tail = FALSE)
    p_wall <- if (x == 0) {
      1
    } else {
      BiasedUrn::pWNCHypergeo(x - 1, n_cat, N - n_cat, n_focal, odds,
                              lower.tail = FALSE)
    }
    list(n_cat = n_cat, x = x, odds = odds, p_hyper = p_hyper, p_wall = p_wall,
         gene_ids = paste(hits, collapse = "/"))
  }
  stats_list <- lapply(sets, per_cat)

  res <- tibble::tibble(
    term = names(sets),
    n_focal_in_cat = vapply(stats_list, \(s) s$x, integer(1)),
    n_cat = vapply(stats_list, \(s) s$n_cat, integer(1)),
    odds = vapply(stats_list, \(s) s$odds, numeric(1)),
    p_uncorrected = vapply(stats_list, \(s) s$p_hyper, numeric(1)),
    p_wall = vapply(stats_list, \(s) s$p_wall, numeric(1)),
    gene_ids = vapply(stats_list, \(s) s$gene_ids, character(1))
  )
  # clusterProfiler::enricher() -- and so run_ora() -- only returns categories
  # containing at least one focal gene. Categories with none have p = 1 and can
  # never be significant, but they do enlarge the BH denominator. Drop them by
  # default so padj is directly comparable with run_ora(). (The goseq package
  # itself keeps them; set drop_empty_cats = FALSE to match that instead.)
  if (drop_empty_cats) res <- dplyr::filter(res, n_focal_in_cat > 0)
  if (nrow(res) == 0) return(NULL)

  p_used <- if (method == "Wallenius") res$p_wall else res$p_uncorrected
  res <- res |>
    dplyr::mutate(
      contrast = fcontrast,
      DE_direction = DE_direction,
      padj = stats::p.adjust(p_used, method = "BH"),
      n_focal = n_focal,
      n_total = N,
      fold_enrich = (n_focal_in_cat / n_focal) / (n_cat / n_total),
      method = method
    )

  # -- Descriptions, ontology, redundancy -------------------------------------
  res <- term_map |>
    dplyr::select(term, dplyr::any_of(c("description", "ontology"))) |>
    dplyr::distinct(term, .keep_all = TRUE) |>
    dplyr::right_join(res, by = "term")

  is_sig <- res$padj < p_enrich & res$n_focal_in_cat >= min_DE_in_cat
  res$redundant <- dplyr::if_else(is_sig, FALSE, NA)
  if (simplify_terms && any(is_sig, na.rm = TRUE)) {
    res$redundant <- flag_redundant_goseq(res, is_sig, simplify_cutoff)
  }

  # -- Mean/median LFC of the focal genes in each category --------------------
  if (!is.null(df) && "lfc" %in% colnames(df)) {
    w_lfc <- res |>
      dplyr::select(term, contrast, DE_direction, gene_ids) |>
      tidyr::separate_longer_delim(cols = gene_ids, delim = "/") |>
      dplyr::left_join(
        dplyr::select(df, gene, contrast, lfc),
        by = dplyr::join_by("gene_ids" == "gene", "contrast"),
        relationship = "many-to-many"
      ) |>
      dplyr::summarize(
        mean_lfc = mean(lfc, na.rm = TRUE),
        median_lfc = stats::median(lfc, na.rm = TRUE),
        .by = c("term", "contrast", "DE_direction")
      )
    res <- dplyr::left_join(res, w_lfc, by = c("term", "contrast", "DE_direction"))
  } else {
    res <- dplyr::mutate(res, mean_lfc = NA_real_, median_lfc = NA_real_)
  }

  if (verbose) {
    cat("Contrast: ", fcontrast, " // DE direction: ", DE_direction,
        " // DEGs (w/ term): ", length(focal_genes), " (", n_focal, ")",
        " // background genes w/ term: ", N,
        " // enriched terms: ", sum(is_sig, na.rm = TRUE),
        if (simplify_terms) {
          paste0(" (", sum(is_sig & !res$redundant, na.rm = TRUE), " non-redundant)")
        } else {
          ""
        },
        "\n", sep = "")
  }

  # run_ora()'s columns first and in its order, goseq extras appended
  res |>
    dplyr::arrange(padj) |>
    dplyr::select(
      term, contrast, DE_direction, n_focal_in_cat, padj, redundant,
      description, ontology, gene_ids, mean_lfc, median_lfc,
      n_focal, n_cat, n_total, fold_enrich,
      p_uncorrected, odds, method
    )
}

# Flag significant terms that simplify() would drop as too similar to a more
# significant term (run_goseq() helper). Mirrors flag_redundant_termmap() in
# ora.R, but works on run_goseq()'s own plain results tibble (columns
# term/description/padj/ontology, and a caller-supplied 'is_sig' that already
# folds in min_DE_in_cat) rather than an enrichResult S4 object thresholded on
# p_enrich alone -- kept separate rather than shared for that reason.
# Prints nothing: run_goseq()'s own progress line already reports both the
# significant and the non-redundant term count ("enriched terms: 6 (5
# non-redundant)"), unlike run_ora(), which builds its line in pieces.
flag_redundant_goseq <- function(res, is_sig, simplify_cutoff) {
  redundant <- dplyr::if_else(is_sig, FALSE, NA)
  simplify_fn <- tryCatch(
    utils::getFromNamespace("simplify_internal", "clusterProfiler"),
    error = function(e) NULL
  )
  if (is.null(simplify_fn)) {
    warning("clusterProfiler:::simplify_internal() not found -- skipping simplify")
    return(redundant)
  }
  sig_rows <- res[which(is_sig), ] |>
    dplyr::transmute(ID = term, Description = description, p.adjust = padj, ontology)
  onts <- unique(sig_rows$ontology[!is.na(sig_rows$ontology)])
  keep <- unlist(lapply(onts, function(ont) {
    rows <- sig_rows |> dplyr::filter(ontology == ont)
    if (nrow(rows) == 0) return(NULL)
    sem <- godata_cached(ont)
    simplify_fn(rows, cutoff = simplify_cutoff, measure = "Wang",
                ontology = ont, semData = sem)$ID
  }))
  redundant[is_sig] <- !res$term[is_sig] %in% keep
  redundant
}

# GOSemSim::godata() is rebuilt from scratch on every call, which is wasteful
# when run_goseq() is mapped over many contrasts and directions. Cache per
# ontology for the lifetime of the session.
.godata_cache <- new.env(parent = emptyenv())
godata_cached <- function(ont) {
  if (is.null(.godata_cache[[ont]])) {
    .godata_cache[[ont]] <- GOSemSim::godata(ont = ont, computeIC = FALSE)
  }
  .godata_cache[[ont]]
}
