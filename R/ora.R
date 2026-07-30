# Fold significance and redundancy into a single 'redundant' column:
#   NA    -- not significant (didn't pass p_enrich/q_enrich/min_DE_in_cat)
#   FALSE -- significant, and kept by simplify() (or simplify_terms == FALSE)
#   TRUE  -- significant, but removed by simplify() as too similar to another term
# Treats NA qvalue (e.g. a failed qvalue::qvalue() call, common with few terms or a
# narrow p-value range) as 'not significant' rather than letting it propagate to NA and
# silently emptying downstream filter(redundant == FALSE) calls in cdotplot().
# 'is_redundant' is the raw TRUE/FALSE flag set by flag_redundant_go()/flag_redundant_termmap()
# (always FALSE when simplify_terms == FALSE), and only meaningful for significant terms.
compute_redundant_col <- function(
  padj,
  qvalue,
  count,
  p_enrich,
  q_enrich,
  min_DE_in_cat,
  is_redundant = FALSE
) {
  significant <- dplyr::coalesce(
    padj < p_enrich & qvalue < q_enrich & count >= min_DE_in_cat,
    FALSE
  )
  dplyr::if_else(significant, is_redundant, NA)
}

# Report the number of significant terms (run_ora() progress message). 'redundant_col' is
# the final 'redundant' column (NA = not significant, FALSE/TRUE = significant, kept/removed
# by simplify()). When simplify_terms == TRUE, also reports the count before redundancy
# filtering (i.e. all significant terms, kept or removed) so the user can see both.
report_sig_counts <- function(redundant_col, simplify_terms) {
  n_sig_total <- sum(!is.na(redundant_col))
  if (simplify_terms) {
    n_sig_kept <- sum(!is.na(redundant_col) & !redundant_col)
    cat(
      " // significant terms:",
      n_sig_kept,
      "(",
      n_sig_total,
      "before redundancy filtering )\n"
    )
  } else {
    cat(" // significant terms:", n_sig_total, "\n")
  }
}

# Warn once if an entire qvalue column is NA (see compute_redundant_col() note above).
warn_if_qvalue_all_na <- function(qvalue) {
  if (length(qvalue) > 0 && all(is.na(qvalue))) {
    message(
      "Note: q-values are all NA (qvalue::qvalue() likely failed) -- ",
      "significance is computed from p.adjust only"
    )
  }
}

# Flag redundant GO terms in an enrichGO() result using clusterProfiler::simplify().
# Adds a logical 'redundant' column to res@result: TRUE for terms that were significant
# (p.adjust < p_enrich) but removed by simplify() as too similar to another, more
# significant term. Every tested term is retained; terms that were never significant get
# redundant = FALSE (simplify never considered them). simplify()/simplify_internal()
# computes an all-vs-all similarity matrix over every row it is given, so we only run it
# on the (typically much smaller) significant subset rather than the full tested set.
# 'res' must be an enrichResult with @ontology already set.
flag_redundant_go <- function(res, p_enrich, simplify_cutoff) {
  res@result$redundant <- FALSE
  full_result <- res@result
  # Deliberately only thresholds on p_enrich (not q_enrich/min_DE_in_cat, which also feed
  # the final 'sig' flag in cp_to_df()/run_ora()) -- simplify() is given a superset of the
  # eventually-'sig' terms, which is fine since redundant is only meaningful for 'sig' rows.
  sig_result <- full_result |> dplyr::filter(p.adjust < p_enrich)
  if (nrow(sig_result) > 0) {
    res@result <- sig_result
    simplified <- clusterProfiler::simplify(
      res,
      cutoff = simplify_cutoff,
      measure = "Wang"
    )
    removed_ids <- setdiff(sig_result$ID, simplified@result$ID)
    full_result$redundant <- full_result$ID %in% removed_ids
  }
  res@result <- full_result
  return(res)
}

# Flag redundant terms in an enricher() (manual term_map) result. Unlike the GO/OrgDb
# path, a term_map result can mix ontologies, so significant terms are grouped by the
# term's ontology and clusterProfiler:::simplify_internal() is run per ontology directly
# on the rows. As with flag_redundant_go(), every tested term is retained and redundant
# significant terms are flagged (not dropped). Returns res with res@result$redundant set.
flag_redundant_termmap <- function(
  res,
  term_map,
  p_enrich,
  simplify_cutoff,
  verbose = TRUE
) {
  res@result$redundant <- FALSE
  ont_lookup <- term_map |>
    dplyr::select(term, ontology) |>
    dplyr::distinct(term, .keep_all = TRUE)
  # Deliberately only thresholds on p_enrich (see the matching note in flag_redundant_go())
  res_df <- res@result |>
    dplyr::left_join(ont_lookup, by = c("ID" = "term")) |>
    dplyr::filter(p.adjust < p_enrich)
  ontologies <- unique(res_df$ontology[!is.na(res_df$ontology)])

  # clusterProfiler:::simplify_internal() is not exported and may change/disappear across
  # clusterProfiler versions -- look it up once and fall back to 'no simplify' if missing,
  # rather than erroring out.
  simplify_fn <- tryCatch(
    utils::getFromNamespace("simplify_internal", "clusterProfiler"),
    error = function(e) NULL
  )
  if (is.null(simplify_fn)) {
    warning(
      "clusterProfiler:::simplify_internal() not found (clusterProfiler version may ",
      "have changed) -- skipping simplify, no terms flagged as redundant"
    )
    return(res)
  }

  simplified_rows <- lapply(ontologies, function(ont) {
    rows <- res_df |> dplyr::filter(ontology == ont)
    if (nrow(rows) == 0) {
      return(NULL)
    }
    semData <- GOSemSim::godata(ont = ont, computeIC = FALSE)
    simplify_fn(
      rows,
      cutoff = simplify_cutoff,
      measure = "Wang",
      ontology = ont,
      semData = semData
    )
  })
  keep_ids <- do.call(rbind, simplified_rows)$ID
  removed_ids <- setdiff(res_df$ID, keep_ids)
  if (verbose) {
    cat(" // simplified:", nrow(res_df), "->", length(keep_ids), "terms")
  }
  res@result$redundant <- res@result$ID %in% removed_ids
  return(res)
}

#' Convert a clusterProfiler enrichment result to a tidy dataframe
#'
#' Converts a `clusterProfiler` enrichment result to a tidy tibble. This is
#' the same conversion [run_ora()] applies internally when `return_df =
#' TRUE`, exposed separately so it can be used on a `clusterProfiler` result
#' produced outside [run_ora()] (e.g. a direct `enricher()`/`enrichGO()`
#' call) and fed to [cdotplot()].
#'
#' @param res A `clusterProfiler` `enrichResult` (S4 object), or a tibble
#'   already built via `purrr::map_dfr()`.
#' @param df Original DE results (needs `gene` + `lfc` columns); when
#'   supplied, adds `mean_lfc`/`median_lfc` per term.
#' @param term_map Optional term-to-gene map; when it has an `ontology`
#'   column, adds an `ontology` column (term_map path).
#' @param contrast Value written into the `contrast` column.
#' @param DE_direction Value written into the `DE_direction` column.
#' @param p_enrich Adjusted p-value threshold for significance (`redundant`
#'   column).
#' @param q_enrich Q-value threshold for significance (`redundant` column).
#' @param min_DE_in_cat Min. number of DE genes in a category for
#'   significance (`redundant` column).
#'
#' @return A tibble with one row per tested term.
#' @export
cp_to_df <- function(
  res,
  df = NULL,
  term_map = NULL,
  contrast = NA,
  DE_direction = NA,
  p_enrich = 0.05,
  q_enrich = 0.2,
  min_DE_in_cat = 2
) {
  fcontrast <- contrast

  # NOTE: when simplify_terms == TRUE upstream, a 'redundant' column will already be present,
  # flagging significant terms simplify() removed as too similar to another, more sig. term.
  # 'res' may already be a tibble here (OrgDb/GO path, built via map_dfr()) or still be
  # the S4 enrichResult object (term_map path) -- handle both.
  if (methods::is(res, "enrichResult")) {
    if (!"redundant" %in% colnames(res@result)) {
      res@result$redundant <- FALSE
    }
  } else if (!"redundant" %in% colnames(res)) {
    res$redundant <- FALSE
  }
  res <- tibble::as_tibble(res)
  warn_if_qvalue_all_na(res$qvalue)
  res <- res |>
    dplyr::mutate(
      redundant = compute_redundant_col(
        p.adjust,
        qvalue,
        Count,
        p_enrich,
        q_enrich,
        min_DE_in_cat,
        redundant
      ),
      contrast = fcontrast,
      DE_direction = DE_direction
    ) |>
    dplyr::select(
      contrast,
      DE_direction,
      term = ID,
      n_focal_in_cat = Count,
      GeneRatio,
      BgRatio,
      padj = p.adjust,
      redundant,
      description = Description,
      dplyr::any_of("ontology"),
      gene_ids = geneID
    )

  # Add ontology info for GO
  if (!is.null(term_map) && "ontology" %in% colnames(term_map)) {
    res <- term_map |>
      dplyr::select(term, ontology) |>
      dplyr::distinct(term, .keep_all = TRUE) |>
      dplyr::right_join(res, by = "term") |>
      dplyr::relocate(ontology, .before = "gene_ids")
  }

  # Add mean & median LFC value
  # NOTE: joined on gene *and* contrast -- 'df' is the full (multi-contrast) DE results
  # table, so a gene-only join would average LFCs across every contrast it appears in.
  if (!is.null(df)) {
    w_lfc <- res |>
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
    res <- dplyr::left_join(
      res,
      w_lfc,
      by = c("term", "contrast", "DE_direction")
    )
  }

  # Add gene numbers and enrichment ratio
  if (nrow(res) == 0) {
    # With 0 rows, separate_wider_delim() can't infer a split width from
    # GeneRatio/BgRatio and would drop them, breaking the mutate() below --
    # so just add the expected (empty, correctly-typed) columns directly
    res <- res |>
      dplyr::select(-GeneRatio, -BgRatio) |>
      dplyr::mutate(
        n_focal = integer(0),
        n_cat = integer(0),
        n_total = integer(0),
        fold_enrich = numeric(0)
      )
  } else {
    res <- res |>
      tidyr::separate_wider_delim(
        cols = c("GeneRatio", "BgRatio"),
        delim = "/",
        names_sep = "_"
      ) |>
      dplyr::mutate(
        # Total nr of DE genes (in+not in the focal term)
        n_focal = as.integer(GeneRatio_2),
        # Total nr of genes in the functional term
        n_cat = as.integer(BgRatio_1),
        # Total nr of genes tested
        n_total = as.integer(BgRatio_2)
      ) |>
      dplyr::select(-GeneRatio_1, -GeneRatio_2, -BgRatio_1, -BgRatio_2) |>
      dplyr::mutate(
        fold_enrich = (n_focal_in_cat / n_focal) / (n_cat / n_total)
      )
  }

  return(res)
}

# Build the focal gene list and background universe for an ORA (run_ora() helper).
# Handles the two input modes: (a) a DE results 'df' (+ contrast/DE_direction), from which
# the DEGs are derived, or (b) a pre-built 'focal_genes' vector. Also builds the background
# 'universe', checks for duplicated gene IDs, prints the progress report, and returns
# everything run_ora() needs downstream. Returns a list with elements: focal_genes,
# univ_vec, DE_direction, skip (logical -- TRUE if enrichment should be skipped).
prep_ora_genes <- function(
  df,
  init_df,
  focal_genes,
  contrast,
  DE_direction,
  term_map,
  universe,
  exclude_nontested,
  allow_dups,
  p_DE = 0.05,
  verbose = TRUE
) {
  fcontrast <- contrast

  if (is.null(focal_genes)) {
    # NOTE: 'df' has already had log2FoldChange -> lfc renamed by run_ora(), up front.
    if (!"isDE" %in% colnames(df)) {
      df <- df |> dplyr::mutate(isDE = ifelse(padj < p_DE, TRUE, FALSE))
    }

    # Filter DE results to only get those for the focal contrast
    df <- df |> dplyr::filter(contrast == fcontrast)
    stopifnot(
      "ERROR: no rows in DE results dataframe after contrast filtering" = nrow(
        df
      ) >
        0
    )

    # Filter the DE results, if needed: only take over- or under-expressed
    if (DE_direction == "up") {
      df <- df |> dplyr::filter(lfc > 0)
    }
    if (DE_direction == "down") {
      df <- df |> dplyr::filter(lfc < 0)
    }

    # Create a vector with DEGs
    focal_genes <- df |> dplyr::filter(isDE) |> dplyr::pull(gene)
  } else {
    DE_direction <- NA_character_
  }

  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(focal_genes))) {
    if (allow_dups) {
      focal_genes <- unique(focal_genes)
    } else {
      stop(
        "Found duplicated gene IDs: you probably have multiple 'contrasts' in your input df"
      )
    }
  }

  # Build the background 'universe' of genes.
  # With a term_map: genes that were tested for DE *and* occur in the term map
  # (Excluding the latter is equivalent to goseq's 'use_genes_without_cat=FALSE',
  # and this is done by default by ClusterProfiler -- but non-tested genes *are* included)
  if (!is.null(term_map)) {
    if (!is.null(universe)) {
      univ_vec <- universe[universe %in% term_map$gene]
    } else if (exclude_nontested == TRUE && !is.null(init_df)) {
      univ_vec <- init_df |>
        dplyr::filter(!is.na(padj), gene %in% term_map$gene) |>
        dplyr::pull(gene) |>
        unique()
    } else {
      univ_vec <- NULL
    }
  } else {
    if (!is.null(universe)) {
      univ_vec <- universe
    } else if (exclude_nontested == TRUE && !is.null(init_df)) {
      univ_vec <- init_df |>
        dplyr::filter(!is.na(padj)) |>
        dplyr::pull(gene) |>
        unique()
    } else {
      univ_vec <- NULL
    }
  }

  # Report progress (mode-specific), and flag whether enrichment should be skipped
  skip <- FALSE
  if (!is.null(term_map)) {
    genes_in_map <- focal_genes[focal_genes %in% term_map$gene]
    # NULL univ_vec means 'use clusterProfiler's default universe', not 'zero background genes'
    univ_label <- if (is.null(univ_vec)) "default" else length(univ_vec)
    if (verbose) {
      cat(
        "Contrast: ",
        fcontrast,
        " // DE direction: ",
        DE_direction,
        " // DEGs (w/ term): ",
        length(focal_genes),
        " (",
        length(genes_in_map),
        ")",
        " // background genes w/ term: ",
        univ_label,
        sep = ""
      )
    }
    if (length(focal_genes) <= 2) {
      message(" // Skipping enrichment analysis: fewer than 3 DEGs\n")
      skip <- TRUE
    } else if (length(genes_in_map) == 0) {
      message(
        "WARNING: None of the DE genes are in the GO/KEGG term dataframe ('term_map')"
      )
      cat("First gene IDs from DE results: ", utils::head(focal_genes), "\n")
      cat("First gene IDs from term_map: ", utils::head(term_map$gene), "\n")
      cat("Skipping enrichment analysis...\n")
      skip <- TRUE
    }
  } else {
    if (verbose) {
      cat(
        "Contrast: ",
        fcontrast,
        " // DE direction: ",
        DE_direction,
        " // Nr DE genes: ",
        length(focal_genes),
        sep = ""
      )
    }
    if (length(focal_genes) <= 1) {
      message("// Skipping enrichment analysis: 1 or 0 DEGs\n")
      skip <- TRUE
    }
  }

  list(
    focal_genes = focal_genes,
    univ_vec = univ_vec,
    DE_direction = DE_direction,
    skip = skip
  )
}

# Run enrichment with a manual term_map via clusterProfiler::enricher() (run_ora() engine).
# Prepares the TERM2GENE / TERM2NAME mappings from term_map, runs enricher(), and (when
# simplify_terms) flags redundant terms. Returns the enrichResult (S4) object, or NULL.
run_enricher <- function(
  focal_genes,
  term_map,
  univ_vec,
  min_cat_size,
  max_cat_size,
  filter_no_descrip,
  simplify_terms,
  p_enrich,
  simplify_cutoff,
  verbose = TRUE
) {
  # Prep term mappings - term-to-gene
  term2gene <- term_map |> dplyr::select(term, gene)

  # Prep term mappings - term-to-name (description)
  if ("description" %in% colnames(term_map)) {
    term2name <- term_map |> dplyr::select(term, description)
    if (filter_no_descrip == TRUE) {
      n_before <- length(unique(term2name$term))
      term2name <- term2name[!is.na(term2name$description), ]
      n_removed <- n_before - length(unique(term2name$term))
      if (n_removed > 0) {
        message("Note: removed ", n_removed, " terms with no description")
      }
      # Also drop those terms from term2gene, or they'd still be tested/returned
      # (with an NA description) despite the message above claiming removal
      term2gene <- term2gene |> dplyr::filter(term %in% term2name$term)
    }
  } else {
    term2name <- NA
  }

  res <- clusterProfiler::enricher(
    gene = focal_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    universe = univ_vec,
    minGSSize = min_cat_size,
    maxGSSize = max_cat_size,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  if (simplify_terms == TRUE && !is.null(res) && nrow(res) > 0) {
    # Always return every tested term; flag (rather than drop) those removed by simplify()
    res@result$redundant <- FALSE
    if ("ontology" %in% colnames(term_map)) {
      res <- flag_redundant_termmap(
        res,
        term_map,
        p_enrich,
        simplify_cutoff,
        verbose = verbose
      )
    } else {
      warning(
        "simplify_terms = TRUE requires an 'ontology' column in term_map -- ",
        "skipping simplify (no terms flagged as redundant)"
      )
    }
  }
  return(res)
}

# Translate/validate a keyType for KEGG functions (shared by run_ora()'s run_enrichkegg()
# and run_gsea()'s KEGG path). KEGG's REST API (queried via KEGG_convert() inside
# enrichKEGG()/gseKEGG()) only recognizes a small, fixed set of gene ID types -- notably
# NOT 'ENSEMBL' -- so an unsupported keyType is caught here with a clear message rather
# than surfacing as a cryptic KEGG HTTP-400/"not supported" error from deep inside
# clusterProfiler. keyType follows enrichGO()'s/gseGO()'s OrgDb-style default
# ('ENTREZID'), translated to KEGG's own vocabulary ('ncbi-geneid') since OrgDb and KEGG
# functions use different keyType names; any other keyType (e.g. 'kegg', 'uniprot') is
# passed through as-is (after validation).
resolve_kegg_keytype <- function(keyType) {
  kegg_keyTypes <- c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot")
  kegg_keyType <- if (identical(keyType, "ENTREZID")) "ncbi-geneid" else keyType
  if (!kegg_keyType %in% kegg_keyTypes) {
    stop(
      "keyType = '",
      keyType,
      "' is not supported by KEGG -- KEGG's API only recognizes: ",
      paste(kegg_keyTypes, collapse = ", "),
      " (or 'ENTREZID', auto-translated to 'ncbi-geneid'). Convert your gene IDs to one ",
      "of these first, e.g. via AnnotationDbi::mapIds(OrgDb, keys, keytype = 'ENSEMBL', ",
      "column = 'ENTREZID')."
    )
  }
  kegg_keyType
}

# Run KEGG enrichment for a model organism via clusterProfiler::enrichKEGG() (run_ora()
# engine), as an alternative to the term_map path. KEGG terms have no GO-like DAG structure,
# so clusterProfiler::simplify() doesn't apply -- simplify_terms is only warned about here,
# never honored; 'redundant' is left FALSE for every term.
run_enrichkegg <- function(
  focal_genes,
  organism,
  keyType,
  univ_vec,
  min_cat_size,
  max_cat_size,
  simplify_terms
) {
  if (simplify_terms == TRUE) {
    warning(
      "simplify_terms = TRUE is not supported for KEGG (no GO-like term structure) -- ",
      "ignoring simplify_terms"
    )
  }
  # Resolved/validated eagerly (as its own statement, not inline as an enrichKEGG() argument):
  # enrichKEGG() only forces its 'keyType' argument deep inside download_KEGG(), after it has
  # already made 1-2 KEGG API calls -- passed inline, an invalid keyType would only surface
  # after those network round-trips instead of failing immediately.
  kegg_keyType <- resolve_kegg_keytype(keyType)
  clusterProfiler::enrichKEGG(
    gene = focal_genes,
    organism = organism,
    keyType = kegg_keyType,
    universe = univ_vec,
    minGSSize = min_cat_size,
    maxGSSize = max_cat_size,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}

# Run GO enrichment with an OrgDb via clusterProfiler::enrichGO() (run_ora() engine).
# When return_df == FALSE, runs a single ontology and returns the enrichResult (S4) object.
# When return_df == TRUE, runs all three ontologies (BP, MF, CC) and returns a combined
# tibble with an 'ontology' column. Flags redundant terms when simplify_terms == TRUE.
run_enrichgo <- function(
  focal_genes,
  OrgDb,
  keyType,
  univ_vec,
  GO_ontology,
  GO_ontologies,
  min_cat_size,
  max_cat_size,
  simplify_terms,
  p_enrich,
  simplify_cutoff,
  return_df
) {
  enrichfun <- function(GO_ontology) {
    clusterProfiler::enrichGO(
      gene = focal_genes,
      OrgDb = OrgDb,
      keyType = keyType,
      universe = univ_vec,
      ont = GO_ontology,
      minGSSize = min_cat_size,
      maxGSSize = max_cat_size,
      pAdjustMethod = "BH",
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }

  if (return_df == FALSE) {
    # If keeping the ClusterProfiler format, can't combine multiple results
    if (is.null(GO_ontology)) {
      GO_ontology <- "BP"
    }
    if (!GO_ontology %in% GO_ontologies) {
      stop(
        "'GO_ontology' must be one of: ",
        paste(GO_ontologies, collapse = ", ")
      )
    }
    warning(
      "return_df = FALSE only runs GO enrichment for a single ontology ('",
      GO_ontology,
      "') -- set return_df = TRUE to run and combine all of: ",
      paste(GO_ontologies, collapse = ", ")
    )
    res <- enrichfun(GO_ontology = GO_ontology)
    # enrichGO() can return NULL (e.g. no term passes minGSSize/maxGSSize)
    if (simplify_terms == TRUE && !is.null(res)) {
      res@ontology <- GO_ontology
      res <- flag_redundant_go(res, p_enrich, simplify_cutoff)
    }
  } else if (simplify_terms == FALSE) {
    # If converting to a df, iterate over the GO ontologies
    res <- purrr::map_dfr(GO_ontologies, function(x) {
      rr <- enrichfun(x)
      if (is.null(rr)) {
        return(NULL)
      }
      tibble::as_tibble(rr) |> dplyr::mutate(ontology = x)
    })
  } else {
    # Simplify each ontology's enrichResult (S4 object) before converting to a tibble.
    res <- purrr::map_dfr(GO_ontologies, function(x) {
      rr <- enrichfun(x)
      if (is.null(rr)) {
        return(NULL)
      }
      rr@ontology <- x
      rr <- flag_redundant_go(rr, p_enrich, simplify_cutoff)
      tibble::as_tibble(rr) |> dplyr::mutate(ontology = x)
    })
  }
  return(res)
}

# Convert a raw DESeq2 results() object (class 'DESeqResults') to the plain DE-results
# tibble run_ora() expects: gene IDs are pulled from rownames, log2FoldChange/baseMean/
# pvalue are renamed to lfc/mean/p, lfcSE and stat are dropped, and rows are sorted by
# padj. A DESeqResults object represents a single comparison and has no 'contrast' column,
# so one is added here (tagged with the value passed to run_ora()'s 'contrast' argument)
# so the existing contrast-filtering logic in prep_ora_genes() applies unchanged.
prep_deseq_df <- function(df, contrast) {
  fcontrast <- if (is.null(contrast)) NA_character_ else contrast
  df |>
    as.data.frame() |>
    tibble::rownames_to_column("gene") |>
    dplyr::rename(
      lfc = log2FoldChange,
      mean = baseMean,
      p = pvalue
    ) |>
    dplyr::select(-dplyr::any_of(c("lfcSE", "stat"))) |>
    dplyr::arrange(padj) |>
    tibble::as_tibble() |>
    dplyr::mutate(contrast = fcontrast)
}

#' Run a GO or KEGG over-representation analysis (ORA)
#'
#' A wrapper around `clusterProfiler::enricher()`/`clusterProfiler::enrichGO()`
#' that derives a focal (DE) gene set and a background universe from a
#' differential expression results table (or a directly supplied gene
#' vector), runs the enrichment analysis, and optionally flags redundant
#' terms and/or converts the result to a tidy tibble.
#'
#' @param df DE results df with at least columns `gene`, `contrast`,
#'   `lfc`/`log2FoldChange`, `padj`. This df should contain *all* genes, not
#'   just significantly DE ones. Can also be a raw DESeq2 `results()` object
#'   (class `DESeqResults`) -- it is converted automatically and tagged with
#'   the `contrast` argument below.
#' @param contrast DE comparison as specified in the `contrast` column in DE
#'   results.
#' @param DE_direction DE direction: `"either"` (all DE genes), `"up"` (lfc >
#'   0), or `"down"` (lfc < 0).
#' @param focal_genes Alternative to `df`: a pre-built vector of focal gene
#'   IDs.
#' @param term_map Manually provide a functional category/term to gene
#'   mapping with columns: 1: term, 2: gene, and optionally 3: description,
#'   4: ontology.
#' @param OrgDb Bioconductor `OrgDb` (alternative to `term_map`, for GO
#'   enrichment in available organisms). Required when `ontology_type =
#'   "GO"`; not used for `ontology_type = "KEGG"` (see `organism`).
#' @param ontology_type `"GO"` (`enrichGO()`, requires `OrgDb`) or `"KEGG"`
#'   (`enrichKEGG()`, requires `organism`). Only applies when using an
#'   `OrgDb`/`organism` instead of a `term_map`.
#' @param organism KEGG species code (e.g. `"hsa"` for human, `"mmu"` for
#'   mouse; see <https://rest.kegg.jp/list/organism> for the full list).
#'   Alternative to `term_map`, used instead of `OrgDb` when `ontology_type =
#'   "KEGG"`.
#' @param keyType OrgDb gene ID type (`ontology_type = "GO"`) or
#'   `enrichKEGG()` `keyType` (`ontology_type = "KEGG"`) -- only applies when
#'   using an `OrgDb`/`organism` instead of a `term_map`. For KEGG, the
#'   default `"ENTREZID"` is translated to `enrichKEGG()`'s own
#'   `"ncbi-geneid"`; any other value is passed through as-is, but must be
#'   one KEGG's own API supports: `"kegg"`, `"ncbi-geneid"`,
#'   `"ncbi-proteinid"`, or `"uniprot"` (notably *not* `"ENSEMBL"` -- convert
#'   Ensembl IDs to one of these first, e.g. via `AnnotationDbi::mapIds()`).
#' @param GO_ontology Only applies when `return_df == FALSE` (default `"BP"`
#'   there, since clusterProfiler can't combine multiple ontologies into a
#'   single `enrichResult`); ignored when `return_df == TRUE`, which always
#'   runs GO analyses across all 3 ontologies (BP, CC, MF).
#' @param p_enrich Adjusted p-value threshold for enrichment.
#' @param q_enrich Q value threshold for enrichment.
#' @param min_DE_in_cat Number-of-DE-genes threshold for enrichment: at least
#'   this number of genes in the ontology category should be DE.
#'   (Occasionally, 'small' categories with 1 DE gene can have p-values below
#'   0.05 -- this excludes those.)
#' @param min_cat_size Min. number of genes in a category (=
#'   clusterProfiler's `minGSSize` argument; note that clusterProfiler's own
#'   default is 10).
#' @param max_cat_size Max. number of genes in a category (=
#'   clusterProfiler's `maxGSSize` argument).
#' @param filter_no_descrip Remove categories/terms with no description (at
#'   least for GO terms, these tend to be old/deprecated ones).
#' @param exclude_nontested Exclude genes that weren't tested (i.e., have
#'   `NA` in the `padj` column) for DE from the `universe` of genes.
#' @param universe Manually specify a background `universe` of genes (vector
#'   with gene IDs).
#' @param allow_dups Allow a gene ID to be present multiple times in a
#'   (single-contrast, single DE-direction) list of DEGs. This should
#'   typically *not* be the case, but could be so when working with gene IDs
#'   (orthologs) from another species than the focal species to run the GO
#'   analysis.
#' @param simplify_terms Whether to use clusterProfiler's `simplify()`
#'   function to flag redundant/similar GO terms among the significant ones.
#'   All tested terms are always returned (regardless of `simplify_terms`)
#'   in a single `redundant` column: `NA` for terms that weren't significant
#'   (`padj < p_enrich`, `qvalue < q_enrich`, `Count >= min_DE_in_cat`),
#'   `FALSE` for significant terms kept by `simplify()` (or when
#'   `simplify_terms = FALSE`), and `TRUE` for significant terms removed by
#'   `simplify()` as too similar to another, more significant term.
#' @param simplify_cutoff Simplify similarity cutoff.
#' @param return_df Convert the result object to a simple dataframe (tibble),
#'   instead of keeping the clusterProfiler object. Should be `FALSE` if you
#'   want to use the enrichPlot functions directly.
#' @param p_DE Adjusted p-value threshold used to derive `isDE` from `df`
#'   when `df` has no `isDE` column.
#' @param verbose Print progress messages (Contrast/DE direction/gene
#'   counts, `simplify()` summary, enriched-term count).
#'
#' @return An `enrichResult` object (`return_df = FALSE`), a tibble
#'   (`return_df = TRUE`), or `NULL` when there were too few (or no mapped)
#'   genes to run the enrichment analysis.
#' @export
run_ora <- function(
  df = NULL,
  contrast = NULL,
  DE_direction = "either",
  focal_genes = NULL,
  term_map = NULL,
  OrgDb = NULL,
  ontology_type = NULL,
  organism = NULL,
  keyType = "ENTREZID",
  GO_ontology = NULL,
  p_enrich = 0.05,
  q_enrich = 0.2,
  min_DE_in_cat = 2,
  min_cat_size = 5,
  max_cat_size = 500,
  filter_no_descrip = TRUE,
  exclude_nontested = TRUE,
  universe = NULL,
  allow_dups = FALSE,
  simplify_terms = FALSE,
  simplify_cutoff = 0.7,
  return_df = FALSE,
  p_DE = 0.05,
  verbose = TRUE
) {
  init_df <- df
  GO_ontologies <- c("BP", "MF", "CC")

  # Validate inputs and determine the analysis mode -----------------------------
  # Exactly one of 'term_map', 'OrgDb', or 'organism' must be supplied
  if (is.null(term_map) && is.null(OrgDb) && is.null(organism)) {
    stop(
      "Supply either a 'term_map', an 'OrgDb' (for GO), or an 'organism' (for KEGG) -- none was given"
    )
  }
  if (!is.null(term_map) && (!is.null(OrgDb) || !is.null(organism))) {
    stop("Supply either a 'term_map' or an 'OrgDb'/'organism', not both")
  }
  mode <- if (!is.null(term_map)) "term_map" else "orgdb"

  # For the OrgDb mode, 'ontology_type' selects the enrichment function
  if (mode == "orgdb") {
    if (is.null(ontology_type) || !ontology_type %in% c("GO", "KEGG")) {
      stop("With an 'OrgDb'/'organism', set ontology_type to 'GO' or 'KEGG'")
    }
    if (ontology_type == "GO" && is.null(OrgDb)) {
      stop("ontology_type = 'GO' requires an 'OrgDb'")
    }
    if (ontology_type == "KEGG" && is.null(organism)) {
      stop(
        "ontology_type = 'KEGG' requires an 'organism' (KEGG species code, e.g. 'hsa')"
      )
    }
  }

  # Need either a DE results 'df' or a pre-built 'focal_genes' vector
  if (is.null(df) && is.null(focal_genes)) {
    stop("Supply either a DE results 'df' or a 'focal_genes' vector")
  }
  if (is.null(focal_genes) && is.null(contrast)) {
    stop("Supply a 'contrast' when providing a DE results 'df'")
  }
  if (!DE_direction %in% c("either", "up", "down")) {
    stop("'DE_direction' must be one of 'either', 'up', or 'down'")
  }

  # Normalize term_map column names once, up front
  if (mode == "term_map") {
    colnames(term_map)[1:2] <- c("term", "gene")
    if (ncol(term_map) > 2) {
      colnames(term_map)[3] <- "description"
    }
    if (ncol(term_map) > 3) {
      colnames(term_map)[4] <- "ontology"
    }
  }

  # A raw DESeq2 results() object needs converting to a plain DE-results tibble before
  # anything downstream can use it (see prep_deseq_df())
  if (!is.null(df) && inherits(df, "DESeqResults")) {
    df <- prep_deseq_df(df, contrast = contrast)
    init_df <- df
  }

  # Normalize the LFC column name once, up front, so every downstream helper (and the
  # 'df' passed into cp_to_df() for mean_lfc/median_lfc) sees a consistent 'lfc' column
  if (
    !is.null(df) &&
      "log2FoldChange" %in% colnames(df) &&
      !"lfc" %in% colnames(df)
  ) {
    df <- df |> dplyr::rename(lfc = log2FoldChange)
    init_df <- df
  }

  # Build the focal gene list and background universe ---------------------------
  genes <- prep_ora_genes(
    df = df,
    init_df = init_df,
    focal_genes = focal_genes,
    contrast = contrast,
    DE_direction = DE_direction,
    term_map = term_map,
    universe = universe,
    exclude_nontested = exclude_nontested,
    allow_dups = allow_dups,
    p_DE = p_DE,
    verbose = verbose
  )
  focal_genes <- genes$focal_genes
  univ_vec <- genes$univ_vec
  DE_direction <- genes$DE_direction
  # NULL only happens when 'focal_genes' was supplied directly with no 'contrast' --
  # keep the eventual 'contrast' column as NA rather than dropping it (mutate(x = NULL) removes it)
  fcontrast <- if (is.null(contrast)) NA_character_ else contrast

  # Skip the enrichment analysis if there are too few (or no mapped) genes
  if (genes$skip) {
    return(NULL)
  }

  # Run the enrichment analysis -------------------------------------------------
  res <- switch(
    mode,
    term_map = run_enricher(
      focal_genes = focal_genes,
      term_map = term_map,
      univ_vec = univ_vec,
      min_cat_size = min_cat_size,
      max_cat_size = max_cat_size,
      filter_no_descrip = filter_no_descrip,
      simplify_terms = simplify_terms,
      p_enrich = p_enrich,
      simplify_cutoff = simplify_cutoff,
      verbose = verbose
    ),
    orgdb = if (ontology_type == "GO") {
      run_enrichgo(
        focal_genes = focal_genes,
        OrgDb = OrgDb,
        keyType = keyType,
        univ_vec = univ_vec,
        GO_ontology = GO_ontology,
        GO_ontologies = GO_ontologies,
        min_cat_size = min_cat_size,
        max_cat_size = max_cat_size,
        simplify_terms = simplify_terms,
        p_enrich = p_enrich,
        simplify_cutoff = simplify_cutoff,
        return_df = return_df
      )
    } else {
      run_enrichkegg(
        focal_genes = focal_genes,
        organism = organism,
        keyType = keyType,
        univ_vec = univ_vec,
        min_cat_size = min_cat_size,
        max_cat_size = max_cat_size,
        simplify_terms = simplify_terms
      )
    }
  )

  # ClusterProfiler may return NULL result for small sets
  if (is.null(res)) {
    message("ClusterProfiler returned NULL (can happen with small sets")
    return(NULL)
  }

  # Process the output
  if (return_df == FALSE) {
    # Every tested term is retained (pvalueCutoff/qvalueCutoff were set to 1 upstream) --
    # fold significance into 'redundant' (NA = not significant, FALSE = significant &
    # kept, TRUE = significant but removed by simplify()), using the same NA-safe rule
    # cp_to_df() uses for the return_df = TRUE path.
    is_redundant <- if ("redundant" %in% colnames(res@result)) {
      res@result$redundant
    } else {
      FALSE
    }
    warn_if_qvalue_all_na(res@result$qvalue)
    res@result$redundant <- compute_redundant_col(
      res@result$p.adjust,
      res@result$qvalue,
      res@result$Count,
      p_enrich,
      q_enrich,
      min_DE_in_cat,
      is_redundant
    )
    if (verbose) {
      report_sig_counts(res@result$redundant, simplify_terms)
    }
  } else {
    # Create a regular df via the standalone conversion helper (see cp_to_df()).
    res <- cp_to_df(
      res,
      df = df,
      term_map = term_map,
      contrast = fcontrast,
      DE_direction = DE_direction,
      p_enrich = p_enrich,
      q_enrich = q_enrich,
      min_DE_in_cat = min_DE_in_cat
    )

    # Report
    if (verbose) {
      report_sig_counts(res$redundant, simplify_terms)
    }
  }

  return(res)
}
