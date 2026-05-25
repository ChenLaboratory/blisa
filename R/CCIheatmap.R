#library(ComplexHeatmap)

#' Heatmap of CCI scores across all ligand-receptor pairs
#'
#' Draws a clustered heatmap (via \code{ComplexHeatmap}) with row annotations
#' for sender and receiver cell types. Row annotations use \code{cell_type_colors},
#' which must be a named character vector of colours defined in the calling
#' environment or passed explicitly.
#'
#' @param CCI_df The \code{CCI_scores} slot of a \code{blisa} object (i.e.
#'   \code{res$CCI_scores}). Must contain columns \code{Sender}, \code{Receiver},
#'   and one column per LR pair.
#' @param include_celltypes Character vector or \code{NULL}. If provided, only
#'   rows where the sender or receiver appears in this vector are kept.
#' @param cell_type_colors Named character vector mapping cell-type names to
#'   colours, used for the sender/receiver row annotations. When \code{NULL}
#'   (default), colours are assigned automatically from the package palette.
#'
#' @return Invisibly returns the \code{Heatmap} object.
#' @export
CCIheatmap <- function(CCI_df, include_celltypes = NULL, cell_type_colors = NULL) {
  # Optional subsetting by cell type
  if (!is.null(include_celltypes)) {
    keep_idx <- CCI_df$Sender %in% include_celltypes | CCI_df$Receiver %in% include_celltypes
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
  lr_cols   <- setdiff(colnames(CCI_df), c("Sender", "Receiver"))
  row_labels <- paste(senders, receivers, sep = " -> ")

  f1 = viridisLite::viridis(10)
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

  # Draw the heatmap with adjusted layout to avoid truncation
  draw(ht,
       heatmap_legend_side = "right",
       padding = unit(c(5, 15, 5, 10), "mm"))

  # Optionally return the heatmap object (useful if chaining)
  invisible(ht)
}
