#' Heatmap of CCI scores across all ligand-receptor pairs
#'
#' Generic function. Draws a clustered heatmap (via \code{ComplexHeatmap})
#' with rows as Sender \eqn{\rightarrow} Receiver cell-type pairs and columns
#' as LR pairs. Row annotations colour-code the sender and receiver cell types.
#'
#' @param x A \code{blisa} object or a CCI scores data frame (the
#'   \code{CCI_scores} slot of a \code{blisa} object). The data frame must
#'   contain columns \code{Sender}, \code{Receiver}, and one column per LR pair.
#' @param ... Additional arguments passed to the relevant method.
#'
#' @return Invisibly returns the \code{Heatmap} object.
#' @seealso \code{\link{plotCCIpair}} for a sender-by-receiver heatmap of a
#'   single LR pair.
#' @export
plotCCI <- function(x, ...) UseMethod("plotCCI")


#' @describeIn plotCCI Method for a \code{blisa} object. Extracts
#'   \code{CCI_scores} and delegates to \code{plotCCI.default}. Stops with an
#'   informative error if \code{CCI_scores} is \code{NULL}.
#'
#' @param include_celltypes Character vector or \code{NULL}. If provided, only
#'   rows where the sender or receiver appears in this vector are kept.
#' @param cell_type_colors Named character vector mapping cell-type names to
#'   colours, used for the sender/receiver row annotations. When \code{NULL}
#'   (default), colours are assigned automatically from the package palette.
#'
#' @export
plotCCI.blisa <- function(x, include_celltypes = NULL,
                          cell_type_colors = NULL, ...) {
  if (is.null(x$CCI_scores))
    stop("CCI_scores is NULL. Run runCCI() first to compute CCI scores.")
  plotCCI.default(x$CCI_scores, include_celltypes = include_celltypes,
                  cell_type_colors = cell_type_colors)
}


#' @describeIn plotCCI Method for a CCI scores data frame (e.g. the
#'   \code{CCI_scores} slot of a \code{blisa} object).
#'
#' @export
plotCCI.default <- function(x, include_celltypes = NULL,
                            cell_type_colors = NULL, ...) {
  CCI_df <- x

  # Optional subsetting by cell type
  if (!is.null(include_celltypes)) {
    keep_idx <- CCI_df$Sender %in% include_celltypes |
                CCI_df$Receiver %in% include_celltypes
    CCI_df <- CCI_df[keep_idx, , drop = FALSE]
    if (nrow(CCI_df) == 0)
      stop("No matching cell type pairs found for the specified include_celltypes.")
  }

  senders   <- CCI_df$Sender
  receivers <- CCI_df$Receiver

  if (is.null(cell_type_colors)) {
    all_cts <- sort(unique(c(senders, receivers)))
    cell_type_colors <- setNames(cols[seq_along(all_cts)], all_cts)
  }

  row_ha <- rowAnnotation(
    Sender = senders,
    Receiver = receivers,
    col = list(
      Sender = cell_type_colors,
      Receiver = cell_type_colors
    ),
    annotation_legend_param = list(Sender = list(show = FALSE),
                                   Receiver = list(show = FALSE)),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side = "top"
  )

  # Score matrix: exclude the Sender / Receiver identifier columns
  lr_cols    <- setdiff(colnames(CCI_df), c("Sender", "Receiver"))
  row_labels <- paste(senders, receivers, sep = " -> ")

  f1 <- viridisLite::viridis(10)
  ht <- ComplexHeatmap::Heatmap(as.matrix(CCI_df[, lr_cols, drop = FALSE]),
                name = "Interaction\nScore",
                col = f1,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                row_labels = row_labels,
                row_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                column_names_rot = 45,
                heatmap_legend_param = list(title_position = "topcenter"),
                left_annotation = row_ha)

  draw(ht,
       heatmap_legend_side = "right",
       padding = unit(c(5, 15, 5, 10), "mm"))

  invisible(ht)
}
