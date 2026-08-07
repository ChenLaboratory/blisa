#' Dot plot ranking LR pairs by number of significant hotspot bins
#'
#' Generic function for ranking LR pairs. Dispatches on the class of \code{x}:
#' \itemize{
#'   \item \code{plotLRrank.blisa} accepts a \code{blisa} object and uses its
#'     \code{LR_results} slot directly.
#'   \item \code{plotLRrank.data.frame} accepts the \code{LR_results} data
#'     frame directly.
#' }
#'
#' @param x A \code{blisa} object or a data frame of LR results. The data frame
#'   must contain columns \code{sig_numbers} and \code{annotation}.
#' @param ... Additional arguments passed to the relevant method.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' # Continuing from the blisa() example:
#' # result <- blisa(spe, bin_size = 50, group = "cell_type")
#' plotLRrank(result, top = 30)
#' plotLRrank(result, top = 20, flip = TRUE)
#' }
#' @export
plotLRrank <- function(x, ...) UseMethod("plotLRrank")


#' @describeIn plotLRrank Method for a \code{blisa} object. Extracts
#'   \code{LR_results} and delegates to \code{plotLRrank.data.frame}.
#'
#' @param top Integer or \code{NULL}. Number of top LR pairs (by
#'   \code{sig_numbers}) to display. Default \code{30}.
#' @param pt_size Numeric. Point size passed to \code{geom_point}. Default 4.
#' @param flip Logical. When \code{TRUE}, LR pairs are placed on the x-axis
#'   and the hotspot count on the y-axis (vertical orientation). Default
#'   \code{FALSE} (LR pairs on y-axis, horizontal orientation).
#'
#' @export
plotLRrank.blisa <- function(x, top = 30, pt_size = 4, flip = FALSE,
                             size_by = "auto", ...) {
  is_pathway <- identical(x$level, "pathway")
  item_label <- if (is_pathway) "Pathway" else "Ligand-Receptor Pair"
  # For pathway results, size points by the number of LR pairs per pathway by
  # default; for LR-level results, no size mapping. Pass size_by explicitly
  # (a column name, or NULL to disable) to override.
  if (identical(size_by, "auto"))
    size_by <- if (is_pathway) "n_LR_pairs" else NULL
  plotLRrank.data.frame(x$LR_results, top = top, pt_size = pt_size, flip = flip,
                        item_label = item_label, size_by = size_by, ...)
}


#' @describeIn plotLRrank Method for a data frame of LR results (e.g. the
#'   \code{LR_results} slot of a \code{blisa} object).
#'
#' @param item_label Character. Axis label for the ranked items. Defaults to
#'   \code{"Ligand-Receptor Pair"}; \code{plotLRrank.blisa} sets it to
#'   \code{"Pathway"} for pathway-level objects.
#' @param size_by Character or \code{NULL}. Name of a numeric column to map to
#'   point size. When \code{NULL}, all points use the fixed \code{pt_size}.
#'   For \code{plotLRrank.blisa} the default is \code{"auto"}, which sizes
#'   pathway-level results (from \code{\link{blisaPathway}}) by
#'   \code{"n_LR_pairs"} and leaves LR-level results unsized; pass a column
#'   name, or \code{NULL} to disable, to override. Relabel the legend with
#'   \code{+ ggplot2::labs(size = ...)}.
#'
#' @export
plotLRrank.data.frame <- function(x, top = 30, pt_size = 4, flip = FALSE,
                                  item_label = "Ligand-Receptor Pair",
                                  size_by = NULL, ...) {
  LR_results <- x

  if (!is.null(top)) {
    LR_results <- LR_results[seq_len(min(top, nrow(LR_results))), ]
  }

  n_shown <- nrow(LR_results)
  title   <- paste0("Top ", n_shown, " by hotspot count")

  LR_results$LR_pair <- rownames(LR_results)
  LR_results <- LR_results[order(-LR_results$sig_numbers), ]
  # flip=FALSE: levels low->high so highest appears at top of y-axis
  # flip=TRUE:  levels high->low so highest appears at left of x-axis
  LR_results$LR_pair <- factor(LR_results$LR_pair,
                            levels = if (flip) LR_results$LR_pair else rev(LR_results$LR_pair))

  # Colour the points by CellChat annotation category when that column is
  # present. Pathway-level tables (from blisaPathway) have no 'annotation'
  # column, so the colour mapping (and legend) are dropped and points are black.
  has_annotation <- "annotation" %in% colnames(LR_results)
  has_size       <- !is.null(size_by) && size_by %in% colnames(LR_results)

  if (has_annotation) {
    # Preset colors for the four standard CellChat categories; any additional
    # annotation values found in the data are assigned colors from the package
    # palette, and NA values fall back to grey.
    known_colors <- c(
      "Secreted Signaling"    = "#90a955",
      "ECM-Receptor"          = "#219ebc",
      "Cell-Cell Contact"     = "#f7b801",
      "Non-protein Signaling" = "#9f86c0"
    )
    extra_levels <- setdiff(
      unique(na.omit(LR_results$annotation)),
      names(known_colors)
    )
    extra_colors <- setNames(cols[seq_along(extra_levels)], extra_levels)
    color_scale  <- scale_color_manual(
      values   = c(known_colors, extra_colors),
      na.value = "grey50"
    )
    color_lab <- "Annotation"
  } else {
    color_scale <- NULL
    color_lab   <- NULL
  }

  # Optionally scale point size by a numeric column (e.g. the number of LR pairs
  # per pathway). When absent, a fixed 'pt_size' is used for all points.
  if (has_size) {
    LR_results$rank_size_val <- LR_results[[size_by]]
    size_scale <- scale_size_continuous(name = size_by)
  } else {
    size_scale <- NULL
  }

  # Build the point aesthetic mapping from the active encodings.
  map <- aes()
  if (has_annotation) map <- modifyList(map, aes(colour = annotation))
  if (has_size)       map <- modifyList(map, aes(size   = rank_size_val))
  point_layer <- if (has_size) geom_point(mapping = map) else
    geom_point(mapping = map, size = pt_size)

  if (flip) {
    p <- ggplot(LR_results, aes(x = LR_pair, y = sig_numbers)) +
      point_layer +
      color_scale +
      size_scale +
      scale_y_continuous(expand = expansion(add = 100)) +
      labs(x = item_label, y = "Sig Spot Numbers",
           color = color_lab, title = title) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title    = element_text(size = 12, face = "bold"),
        legend.text     = element_text(size = 10),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text.x     = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y     = element_text(angle = 0, hjust = 1)
      ) +
      coord_cartesian(clip = "off")
  } else {
    p <- ggplot(LR_results, aes(x = sig_numbers, y = LR_pair)) +
      point_layer +
      color_scale +
      size_scale +
      scale_x_continuous(expand = expansion(add = 100)) +
      labs(x = "Sig Spot Numbers", y = item_label,
           color = color_lab, title = title) +
      theme_minimal() +
      theme(
        legend.position = "right",
        legend.title    = element_text(size = 12, face = "bold"),
        legend.text     = element_text(size = 10),
        panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        axis.text.y     = element_text(size = 12),
        axis.text.x     = element_text(angle = 45, hjust = 1)
      ) +
      coord_cartesian(clip = "off")
  }

  return(p)
}
