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
#' @param top_lr Integer or \code{NULL}. Number of top-ranked LR pairs (by
#'   \code{sig_numbers}) to display as columns. LR pairs in \code{CCI_scores}
#'   are already ordered by rank, so this simply takes the first \code{top_lr}
#'   columns. Default \code{20}.
#' @param top_pairs Integer or \code{NULL}. Number of top sender-receiver pairs
#'   to display as rows, ranked by their maximum interaction score across the
#'   displayed LR pairs (after \code{top_lr} is applied). When \code{NULL} all
#'   rows are shown. Default \code{30}.
#' @param sender Character vector or \code{NULL}. If provided, only rows where
#'   \code{Sender} is in this vector are kept. Applied independently of
#'   \code{receiver} (AND logic when both are supplied). Default \code{NULL}
#'   (all senders).
#' @param receiver Character vector or \code{NULL}. If provided, only rows
#'   where \code{Receiver} is in this vector are kept. Applied independently of
#'   \code{sender} (AND logic when both are supplied). Default \code{NULL}
#'   (all receivers).
#' @param cell_type_colors Named character vector mapping cell-type names to
#'   colours, used for the sender/receiver row annotations. When \code{NULL}
#'   (default), colours are assigned automatically from the package palette.
#'
#' @export
plotCCI.blisa <- function(x, top_lr = 20, top_pairs = 30,
                          sender = NULL, receiver = NULL,
                          cell_type_colors = NULL, ...) {
  if (is.null(x$CCI_scores))
    stop("CCI_scores is NULL. Run runCCI() first to compute CCI scores.")
  plotCCI.default(x$CCI_scores, top_lr = top_lr, top_pairs = top_pairs,
                  sender = sender, receiver = receiver,
                  cell_type_colors = cell_type_colors)
}


#' @describeIn plotCCI Method for a CCI scores data frame (e.g. the
#'   \code{CCI_scores} slot of a \code{blisa} object).
#'
#' @export
plotCCI.default <- function(x, top_lr = 20, top_pairs = 30,
                            sender = NULL, receiver = NULL,
                            cell_type_colors = NULL, ...) {
  CCI_df <- x

  # Filter rows by sender (AND logic with receiver)
  if (!is.null(sender)) {
    CCI_df <- CCI_df[CCI_df$Sender %in% sender, , drop = FALSE]
    if (nrow(CCI_df) == 0)
      stop("No rows remaining after filtering by sender.")
  }

  # Filter rows by receiver (AND logic with sender)
  if (!is.null(receiver)) {
    CCI_df <- CCI_df[CCI_df$Receiver %in% receiver, , drop = FALSE]
    if (nrow(CCI_df) == 0)
      stop("No rows remaining after filtering by receiver.")
  }

  # Subset to top_lr LR pair columns (already ranked by sig_numbers)
  lr_cols <- setdiff(colnames(CCI_df), c("Sender", "Receiver"))
  if (!is.null(top_lr))
    lr_cols <- lr_cols[seq_len(min(top_lr, length(lr_cols)))]

  # Subset to top_pairs rows by max score across the displayed LR columns
  if (!is.null(top_pairs) && nrow(CCI_df) > top_pairs) {
    score_mat <- as.matrix(CCI_df[, lr_cols, drop = FALSE])
    row_max   <- apply(score_mat, 1, max, na.rm = TRUE)
    keep      <- order(row_max, decreasing = TRUE)[seq_len(top_pairs)]
    CCI_df    <- CCI_df[keep, , drop = FALSE]
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
