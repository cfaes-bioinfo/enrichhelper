#' Cleveland dotplot of enrichment results
#'
#' Plots a tidy enrichment-results tibble (e.g. from [run_ora()] or
#' [cp_to_df()] with `return_df = TRUE`) as a Cleveland dotplot, with one dot
#' per row of `df`. `cdotplot()` does not filter for significance or
#' redundancy itself -- filter `df` (e.g. on a `redundant` column: `NA` =
#' not significant, `FALSE` = significant & kept, `TRUE` = significant but
#' removed by `simplify()`) before calling this function.
#'
#' @param df Dataframe with enrichment results from [run_ora()], already
#'   filtered to the rows to plot (e.g. significant, non-redundant terms).
#' @param contrasts One or more contrasts (default: all).
#' @param DE_dirs One or more DE directions (default: all).
#' @param x_var Column in `df` to plot along the x axis (`"padj_log"` will
#'   be computed from `padj`).
#' @param fill_var Column in `df` to vary fill color by (`"padj_log"` will
#'   be computed from `padj`).
#' @param label_var Column in `df` with a number to add as a label in the
#'   circles (optional; when not provided, no labels are shown).
#' @param facet_var1 Column in `df` to facet by.
#' @param facet_var2 Second column in `df` to facet by (e.g. `"ontology"`
#'   for GO).
#' @param facet_to_columns When only using one `facet_var1`, facets are
#'   columns (or rows).
#' @param facet_scales Facet scales: `"fixed"`, `"free"`, `"free_x"`, or
#'   `"free_y"`.
#' @param facet_label_fun Function to label facets with.
#' @param x_title X-axis title.
#' @param ylab_size Size of y-axis labels (= term labels).
#' @param add_term_id Add term ID (e.g., `"GO:0009539"`) to its description.
#' @param point_size Point size.
#' @param label_chars Truncate the term labels to this many characters (+10
#'   when including the term ID).
#' @param scico_palette Name of a `scico::scale_color_scico()` palette
#'   (e.g. `"batlow"`) to use for a numeric `fill_var` instead of the
#'   default viridis scale (ignored for `fill_var %in% c("mean_lfc",
#'   "median_lfc")`, which always use the diverging scale).
#' @param sort_fun Function used to aggregate `abs(x_var)` across a term's
#'   row(s) (a term can have more than one row when faceting) into a single
#'   value that determines its position, with the most extreme term at the
#'   top. Default `max`; use `mean` for an ordering that is less sensitive
#'   to a single extreme facet value. Must accept an `na.rm` argument (as
#'   `max()` and `mean()` do).
#'
#' @return A `ggplot` object.
#' @export
cdotplot <- function(
  df,
  contrasts = NULL,
  DE_dirs = NULL,
  x_var = "padj_log",
  fill_var = "median_lfc",
  label_var = NULL,
  facet_var1 = NULL,
  facet_var2 = NULL,
  facet_to_columns = TRUE,
  facet_scales = NULL,
  facet_label_fun = "label_value",
  x_title = NULL,
  ylab_size = 10,
  add_term_id = FALSE,
  point_size = 6,
  label_chars = 60,
  scico_palette = NULL,
  sort_fun = max
) {
  # Constants
  y_var <- "term"

  # Check that columns needed unconditionally are present, with a clear message
  # naming what's missing, rather than an opaque error deep inside dplyr/ggplot
  required_cols <- c("contrast", "padj", "term", "description")
  missing_required <- setdiff(required_cols, colnames(df))
  if (length(missing_required) > 0) {
    stop(
      "df is missing required column(s): ",
      paste(missing_required, collapse = ", ")
    )
  }
  if (!is.null(DE_dirs) && !"DE_direction" %in% colnames(df)) {
    stop(
      "df is missing the 'DE_direction' column, needed because 'DE_dirs' was specified"
    )
  }
  if (!is.null(facet_var2) && is.null(facet_var1)) {
    stop("facet_var1 must be specified when facet_var2 is used")
  }
  requested_cols <- c(x_var, fill_var, label_var, facet_var1, facet_var2)
  requested_cols <- setdiff(requested_cols, c("padj_log", required_cols))
  missing_requested <- setdiff(requested_cols, colnames(df))
  if (length(missing_requested) > 0) {
    stop(
      "df is missing column(s) requested via x_var/fill_var/label_var/facet_var1/",
      "facet_var2: ",
      paste(missing_requested, collapse = ", ")
    )
  }

  # Select contrasts & DE directions
  if (is.null(contrasts)) {
    contrasts <- unique(df$contrast)
  }

  # Prep the df
  # NOTE: cdotplot() does not filter for significance/redundancy itself -- filter 'df'
  # (e.g. on 'redundant': NA = not significant, FALSE = significant & kept, TRUE =
  # significant but removed by simplify()) before calling this function.
  df <- df |>
    dplyr::filter(contrast %in% contrasts) |>
    dplyr::mutate(padj_log = -log10(padj))
  if (!is.null(DE_dirs)) {
    df <- df |> dplyr::filter(DE_direction %in% DE_dirs)
  }
  if (nrow(df) == 0) {
    stop(
      "No rows left to plot after filtering for the given contrasts/DE_dirs"
    )
  }

  # Order terms by sort_fun(abs(x_var)), most extreme at the top
  # A term can appear in multiple rows (e.g. once per facet), but the term
  # axis is shared across facets, so we need a single, term-level order
  term_order <- df |>
    dplyr::group_by(term) |>
    dplyr::summarise(.sort_key = sort_fun(abs(.data[[x_var]]), na.rm = TRUE)) |>
    dplyr::arrange(.sort_key) |>
    dplyr::pull(term)
  df <- df |> dplyr::mutate(term = factor(term, levels = term_order))

  # Modify the term description
  df <- df |>
    dplyr::mutate(
      # If there is no description, use the term ID (must run before
      # capitalization below -- paste0() turns NA into the string "NA", which
      # would otherwise defeat this is.na() check)
      description = ifelse(is.na(description), term, description),
      # Capitalize the first letter
      description = paste0(
        toupper(substr(description, 1, 1)),
        substr(description, 2, nchar(description))
      )
    )
  if (add_term_id == TRUE) {
    df <- df |> dplyr::mutate(description = paste0(term, " - ", description))
    label_chars <- label_chars + 10
  }
  df <- df |>
    dplyr::mutate(
      description = stringr::str_trunc(description, width = label_chars)
    )

  # Create a label lookup for the term description
  # The problem is that abbreviated terms could be non-unique!
  label_df <- df |>
    dplyr::distinct(term, .keep_all = TRUE) |>
    dplyr::select(term, description)
  label_lookup_vec <- label_df$description
  names(label_lookup_vec) <- label_df$term

  # Legend position and title
  if (x_var == fill_var) {
    legend_pos <- "none"
  } else {
    legend_pos <- "top"
  }
  if (fill_var == "median_lfc") {
    color_name <- expression(paste("Median log"[2] * "-fold change of genes"))
  } else if (fill_var == "mean_lfc") {
    color_name <- expression(paste("Mean log"[2] * "-fold change of genes"))
  } else if (fill_var == "padj_log") {
    color_name <- expression("-Log"[10] * " P")
  } else {
    color_name <- fill_var
  }

  # X-axis title (only fill in a default -- don't clobber a user-supplied x_title)
  if (is.null(x_title)) {
    if (x_var == "padj_log") {
      x_title <- expression("-Log"[10] * " P")
    } else if (x_var == "padj") {
      x_title <- "Adjusted p-value of term"
    } else if (x_var == "fold_enrich") {
      x_title <- "Fold enrichment of term"
    } else if (x_var == "median_lfc") {
      x_title <- expression(paste("Median log"[2] * "-fold change of genes"))
    } else if (x_var == "mean_lfc") {
      x_title <- expression(paste("Mean log"[2] * "-fold change of genes"))
    }
  }

  # Color scale - https://carto.com/carto-colors/
  if (fill_var %in% c("mean_lfc", "median_lfc")) {
    col_scale <- colorspace::scale_color_continuous_divergingx(
      palette = "Tropic",
      mid = 0.0,
      na.value = "grey97",
      name = color_name,
      rev = TRUE
    )
  } else if (is.numeric(df[[fill_var]])) {
    if (!is.null(scico_palette)) {
      col_scale <- scico::scale_color_scico(
        palette = scico_palette,
        na.value = "grey95",
        name = color_name
      )
    } else {
      col_scale <- ggplot2::scale_color_viridis_c(
        option = "D",
        na.value = "grey95",
        name = color_name,
      )
    }
  } else {
    col_scale <- ggplot2::scale_color_brewer(palette = "Dark2")
  }

  # X-axis left-hand expansion
  if (x_var %in% c("padj_log", "fold_enrich")) {
    expand_min <- 0.02
  } else {
    expand_min <- 0.12
  }

  # Create the base plot
  p <- ggpubr::ggdotchart(
    df,
    x = y_var,
    y = x_var,
    label = label_var,
    color = fill_var,
    sorting = "none", # Keep df's row order as-is (no sorting)
    add = "segments", # Add segments from y = 0 to dots
    rotate = TRUE, # Rotate vertically
    dot.size = point_size,
    font.label = list(color = "white", size = point_size + 2, vjust = 0.5),
    ggtheme = ggplot2::theme_bw()
  )

  # Formatting
  p <- p +
    ggplot2::scale_x_discrete(labels = label_lookup_vec) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(expand_min, 0.12))
    ) +
    col_scale +
    ggplot2::labs(x = NULL) +
    ggplot2::theme(
      legend.position = legend_pos,
      plot.margin = ggplot2::margin(0.5, 0.5, 0.5, 0.5, unit = "cm"),
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      strip.text.y = ggplot2::element_text(angle = 270, face = "bold"),
      strip.placement = "outside",
      axis.title.x = ggplot2::element_text(
        size = 12,
        margin = ggplot2::margin(t = 0.5, b = 0.5, unit = "cm")
      ),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 10),
      axis.text.y = ggplot2::element_text(size = ylab_size),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_blank(),
      panel.spacing = grid::unit(11, "pt") # 2x the theme_bw() default
    )

  # Other formatting
  if (!is.null(x_title)) {
    p <- p + ggplot2::labs(y = x_title)
  }

  # Vertical lines at 0 (only if x_var actually has values below 0)
  if (any(df[[x_var]] < 0, na.rm = TRUE)) {
    p <- p +
      ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 1)

    # Find the hline layer and move it to the beginning
    hline_idx <- which(sapply(p$layers, function(x) {
      inherits(x$geom, "GeomHline")
    }))
    if (length(hline_idx) > 0) {
      p$layers <- c(p$layers[hline_idx], p$layers[-hline_idx])
    }
  }

  # Faceting
  if (!is.null(facet_var1)) {
    if (!is.null(facet_var2)) {
      # With 2 faceting variables, use facet_grid()
      if (is.null(facet_scales)) {
        facet_scales <- "free_y"
      }
      p <- p +
        ggplot2::facet_grid(
          rows = ggplot2::vars(.data[[facet_var1]]),
          cols = ggplot2::vars(.data[[facet_var2]]),
          scales = facet_scales,
          space = "free_y",
          labeller = facet_label_fun
        )
    } else if (facet_to_columns) {
      # 1 faceting variable default: facet into columns with facet_row()
      if (is.null(facet_scales)) {
        facet_scales <- "free_x"
      }
      p <- p +
        ggforce::facet_row(
          facets = ggplot2::vars(.data[[facet_var1]]),
          scales = facet_scales,
          space = "free",
          labeller = facet_label_fun
        )
    } else {
      # 1 faceting variable alternative: facet into rows with facet_col()
      if (is.null(facet_scales)) {
        facet_scales <- "free_y"
      }
      p <- p +
        ggforce::facet_col(
          facets = ggplot2::vars(.data[[facet_var1]]),
          scales = facet_scales,
          space = "free",
          labeller = facet_label_fun
        )
    }
  }

  return(p)
}
