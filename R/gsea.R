#' Run a Gene Set Enrichment Analysis (GSEA)
#'
#' A wrapper around `clusterProfiler::GSEA()`/`clusterProfiler::gseGO()`/
#' `clusterProfiler::gseKEGG()` that ranks genes by log2 fold change for a
#' given contrast and runs the GSEA algorithm against either a manually
#' supplied term map or a Bioconductor `OrgDb`/KEGG `organism`.
#'
#' @param df Differential expression results; should at least have columns
#'   `contrast`, `gene`, `lfc`/`log2FoldChange`.
#' @param contrast Contrast ID, should be one of the values in the
#'   `contrast` column in `df`.
#' @param term_map Gene-to-term map (e.g., GO or KEGG). Either `term_map`
#'   (use a 'manual' lookup table), an `OrgDb` (for GO), or an `organism`
#'   (for KEGG) is required.
#' @param OrgDb Bioconductor `OrgDb` (alternative to `term_map`, for GO GSEA
#'   in available organisms). Required when `ontology_type = "GO"`; not used
#'   for `ontology_type = "KEGG"` (see `organism`).
#' @param keyType OrgDb gene ID type (`ontology_type = "GO"`) or
#'   `gseKEGG()` `keyType` (`ontology_type = "KEGG"`) -- only applies when
#'   using an `OrgDb`/`organism` instead of a `term_map`. For KEGG, the
#'   default `"ENTREZID"` is translated to `gseKEGG()`'s own
#'   `"ncbi-geneid"`; any other value is passed through as-is, but must be
#'   one KEGG's own API supports: `"kegg"`, `"ncbi-geneid"`,
#'   `"ncbi-proteinid"`, or `"uniprot"` (notably *not* `"ENSEMBL"` -- convert
#'   Ensembl IDs to one of these first, e.g. via `AnnotationDbi::mapIds()`).
#' @param ontology_type `"GO"` (`gseGO()`, requires `OrgDb`) or `"KEGG"`
#'   (`gseKEGG()`, requires `organism`). Only applies when using an
#'   `OrgDb`/`organism` instead of a `term_map`.
#' @param organism KEGG species code (e.g. `"hsa"` for human, `"mmu"` for
#'   mouse; see <https://rest.kegg.jp/list/organism> for the full list).
#'   Alternative to `term_map`, used instead of `OrgDb` when `ontology_type =
#'   "KEGG"`.
#' @param GO_ontology Only for GO analysis with an `OrgDb` *and* when
#'   `return_df == TRUE` (default is `"BP"` there). In all other cases, GO
#'   analyses are run across all 3 ontologies (BP, CC, MF). (The reason this
#'   applies for GO with an OrgDb while keeping the clusterProfiler format is
#'   that clusterProfiler errors when trying to run `enrichGO()` with all
#'   ontologies at once.)
#' @param p_enrich Adjusted p-value threshold for enrichment.
#' @param simplify_terms Whether to use clusterProfiler's `simplify()`
#'   function to flag redundant/similar GO terms among the significant ones.
#'   With an `OrgDb`, this works directly (GSEA is always run across all 3
#'   GO ontologies). With a `term_map`, this requires an `ontology` column
#'   in `term_map`. Not supported for `ontology_type = "KEGG"` (no GO-like
#'   term structure) -- a warning is emitted and ignored there. Adds a
#'   `redundant` column: `NA` for terms that weren't significant (`padj <
#'   p_enrich`), `FALSE` for significant terms kept by `simplify()` (or when
#'   `simplify_terms = FALSE`), and `TRUE` for significant terms removed by
#'   `simplify()` as too similar to another, more significant term. When
#'   `return_df = FALSE`, terms flagged as redundant are dropped from the
#'   returned object (in addition to non-significant terms).
#' @param simplify_cutoff Simplify similarity cutoff.
#' @param n_perm Number of permutations for GSEA (ClusterProfiler's nPermSimple). Default is 10,000.
#' @param allow_dups Allow a gene ID to be present multiple times in a
#'   (single-contrast) list of DEGs. This should typically *not* be the
#'   case, but could be so when working with gene IDs (orthologs) from
#'   another species than the focal species. NOTE: this computes the mean
#'   LFC across duplicated genes.
#' @param return_df Convert the result object to a simple dataframe
#'   (tibble).
#'
#' @return A GSEA result object (`return_df = FALSE`) or a tibble
#'   (`return_df = TRUE`).
#' @export
run_gsea <- function(
  df,
  contrast,
  term_map = NULL,
  OrgDb = NULL,
  keyType = "ENTREZID",
  ontology_type = NULL,
  organism = NULL,
  GO_ontology = NULL,
  p_enrich = 0.05,
  simplify_terms = FALSE,
  simplify_cutoff = 0.7,
  n_perm = 10000,
  allow_dups = FALSE,
  return_df = FALSE
) {
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
  # For the OrgDb mode, 'ontology_type' selects the enrichment function
  if (is.null(term_map)) {
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

  # Check for name of lfc column, and presence of isDE column
  if ("log2FoldChange" %in% colnames(df) & !"lfc" %in% colnames(df)) {
    df <- df |> dplyr::rename(lfc = log2FoldChange)
  }

  # Prep the df to later create a gene vector
  fcontrast <- contrast
  gene_df <- df |>
    dplyr::filter(contrast == fcontrast, !is.na(lfc)) |>
    dplyr::arrange(dplyr::desc(lfc))
  stopifnot(
    "Error: no rows in DE results dataframe after contrast filtering" = nrow(
      gene_df
    ) >
      0
  )
  n_DE <- sum(gene_df$padj < 0.05, na.rm = TRUE)

  # Check if genes are present multiple times -- this would indicate there are multiple contrasts
  if (any(duplicated(gene_df$gene))) {
    # NOTE: This will compute the mean LFC across duplicated genes!
    if (allow_dups) {
      gene_df <- gene_df |>
        dplyr::summarize(lfc = mean(lfc), .by = gene) |>
        dplyr::arrange(dplyr::desc(lfc))
    } else {
      stop(
        "Found duplicated gene IDs: you probably have multiple 'contrasts' in your input df"
      )
    }
  }

  # Create a vector with lfc's and gene IDs
  lfc_vec <- gene_df$lfc
  names(lfc_vec) <- gene_df$gene

  # Prepare the functional term map and report
  if (!is.null(term_map)) {
    # Rename term_map columns
    colnames(term_map)[1:2] <- c("term", "gene")
    if (ncol(term_map) > 2) {
      colnames(term_map)[3] <- "description"
    }
    if (ncol(term_map) > 3) {
      colnames(term_map)[4] <- "ontology"
    }

    # Prep term mappings - if there's a third column, make a term2name df as well
    term2gene <- term_map[, 1:2]
    term2name <- NA
    if (ncol(term_map) > 2) {
      term2name <- term_map[, c(1, 3)]
    }

    # Check & report
    genes_in_map <- names(lfc_vec)[names(lfc_vec) %in% term2gene[[2]]]
    cat(
      "Contrast: ",
      fcontrast,
      " // Nr genes with term assigned: ",
      length(genes_in_map),
      sep = ""
    )

    if (length(genes_in_map) == 0) {
      message("\nERROR: None of the genes are in the term_map dataframe")
      cat("First gene IDs from DE results: ", utils::head(names(lfc_vec)), "\n")
      cat("First gene IDs from term_map: ", utils::head(term2gene[[2]]), "\n")
      stop()
    }
  } else {
    cat("Contrast:", fcontrast)
  }

  # Set random seed if it doesn't exist
  if (!exists(".Random.seed")) {
    message("Note: no random seed set, setting seed to 1 with `set.seed(1)`")
    set.seed(1)
  }

  # Run the enrichment analysis
  if (!is.null(term_map)) {
    gsea_res <- clusterProfiler::GSEA(
      geneList = lfc_vec,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      nPermSimple = n_perm,
      verbose = FALSE,
      eps = 0,
      seed = TRUE
    )
  } else if (ontology_type == "KEGG") {
    # Resolved/validated eagerly (as its own statement, not inline as a gseKEGG() argument):
    # gseKEGG() only forces its 'keyType' argument deep inside download_KEGG(), after it has
    # already made 1-2 KEGG API calls -- passed inline, an invalid keyType would only surface
    # after those network round-trips instead of failing immediately.
    kegg_keyType <- resolve_kegg_keytype(keyType)
    gsea_res <- clusterProfiler::gseKEGG(
      geneList = lfc_vec,
      organism = organism,
      keyType = kegg_keyType,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      nPermSimple = n_perm,
      verbose = FALSE,
      eps = 0,
      seed = TRUE
    )
  } else {
    gsea_res <- clusterProfiler::gseGO(
      geneList = lfc_vec,
      ont = "ALL",
      OrgDb = OrgDb,
      keyType = keyType,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      nPermSimple = n_perm,
      verbose = FALSE,
      eps = 0,
      seed = TRUE
    )
  }

  # Flag redundant terms (simplify_terms == TRUE), then fold significance and redundancy
  # into a single 'redundant' column: NA = not significant, FALSE = significant & kept
  # (or simplify_terms == FALSE), TRUE = significant but removed by simplify()
  if (simplify_terms) {
    if (!is.null(term_map)) {
      if ("ontology" %in% colnames(term_map)) {
        gsea_res <- flag_redundant_termmap(gsea_res, term_map, p_enrich, simplify_cutoff)
      } else {
        warning(
          "simplify_terms = TRUE requires an 'ontology' column in term_map -- ",
          "skipping simplify (no terms flagged as redundant)"
        )
        gsea_res@result$redundant <- FALSE
      }
    } else if (ontology_type == "KEGG") {
      warning(
        "simplify_terms = TRUE is not supported for KEGG (no GO-like term structure) -- ",
        "ignoring simplify_terms"
      )
      gsea_res@result$redundant <- FALSE
    } else {
      # gseGO(ont = "ALL") sets setType = "GOALL", which simplify() requires
      gsea_res <- flag_redundant_go(gsea_res, p_enrich, simplify_cutoff)
    }
  } else {
    gsea_res@result$redundant <- FALSE
  }
  gsea_res@result$redundant <- dplyr::if_else(
    gsea_res@result$p.adjust < p_enrich,
    gsea_res@result$redundant,
    NA
  )

  # Report
  report_sig_counts(gsea_res@result$redundant, simplify_terms)

  # Return a df, if requested
  if (return_df == FALSE) {
    gsea_res <- gsea_res |> dplyr::filter(p.adjust < p_enrich, !redundant)
  } else {
    # Convert to a df
    gsea_res <- tibble::as_tibble(gsea_res)
    if ("ONTOLOGY" %in% colnames(gsea_res)) {
      gsea_res <- gsea_res |> dplyr::rename(ontology = ONTOLOGY)
    }
    gsea_res <- gsea_res |>
      dplyr::mutate(
        sig = ifelse(p.adjust < p_enrich, TRUE, FALSE),
        contrast = fcontrast,
        # Total nr of genes in the functional term (that were also part of the
        # ranked gene list); total nr of genes tested (i.e. ranked/used for GSEA)
        n_cat = setSize,
        n_total = length(lfc_vec)
      ) |>
      dplyr::select(
        contrast,
        term = ID,
        padj = p.adjust,
        sig,
        redundant,
        description = Description,
        dplyr::any_of("ontology"),
        n_cat,
        n_total,
        gene_ids = core_enrichment
      )

    # Add GO info
    if ("ontology" %in% colnames(term_map)) {
      gsea_res <- term_map |>
        dplyr::select(term, ontology) |>
        dplyr::distinct(term, .keep_all = TRUE) |>
        dplyr::right_join(gsea_res, by = "term") |>
        dplyr::relocate(ontology, .before = "gene_ids")
    }

    # Add mean & median LFC value
    w_lfc <- gsea_res |>
      tidyr::separate_longer_delim(cols = gene_ids, delim = "/") |>
      dplyr::rename(gene = gene_ids) |>
      dplyr::left_join(
        df |> dplyr::select(gene, contrast, lfc),
        by = c("gene", "contrast"),
        relationship = "many-to-many"
      ) |>
      dplyr::summarize(
        mean_lfc = mean(lfc, na.rm = TRUE),
        median_lfc = stats::median(lfc, na.rm = TRUE),
        .by = c("term", "contrast")
      )

    gsea_res <- dplyr::left_join(gsea_res, w_lfc, by = c("term", "contrast"))
  }

  return(gsea_res)
}
