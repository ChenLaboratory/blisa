#library(ComplexHeatmap)

#' plot heatmap of cell-cell interaction scores for all LR pairs
#'
#' @param CCI_df dataframe of cell-cell interaction scores from CCI anlaysis
#' @export
#'
CCIheatmap <- function(CCI_df, include_celltypes = NULL) {
  # Optional subsetting
  if (!is.null(include_celltypes)) {
    # Parse sender and receiver from rownames
    pairs <- strsplit(rownames(CCI_df), "->")
    senders <- sapply(pairs, `[`, 1)
    receivers <- sapply(pairs, `[`, 2)

    # Keep rows where either sender or receiver is in include_celltypes
    keep_idx <- senders %in% include_celltypes | receivers %in% include_celltypes
    CCI_df <- CCI_df[keep_idx, , drop = FALSE]

    if (nrow(CCI_df) == 0) {
      stop("No matching cell type pairs found for the specified include_celltypes.")
    }
  }

  # Row annotations for sender and receiver
  pairs <- strsplit(rownames(CCI_df), "->")
  senders <- sapply(pairs, `[`, 1)
  receivers <- sapply(pairs, `[`, 2)

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

  f1 = viridisLite::viridis(10)
  ht <- ComplexHeatmap::Heatmap(as.matrix(CCI_df),
                name = "Interaction\nScore",
                col = f1,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
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
