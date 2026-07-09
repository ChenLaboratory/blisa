#' Sender-by-receiver heatmap of aggregated CCI scores across LR pairs
#'
#' Generic function. Aggregates CCI scores across all (or the top-ranked)
#' ligand-receptor pairs and draws a clustered receiver-by-sender heatmap, one
#' cell per Sender \eqn{\rightarrow} Receiver combination.
#'
#' @param x A \code{blisa} object or a CCI scores data frame (the
#'   \code{CCI_scores} slot of a \code{blisa} object).
#' @param ... Additional arguments passed to the relevant method.
#'
#' @return A \code{Heatmap} object.
#' @seealso \code{\link{plotCCILR}} for a per-LR-pair version of this plot;
#'   \code{\link{plotCCI}} for a heatmap with LR pairs as columns.
#' @examples
#' \dontrun{
#' # Continuing from the blisa() example:
#' # result <- blisa(spe, bin_size = 50, group = "cell_type")
#' plotCCIsummary(result)
#' plotCCIsummary(result, top_lr = 10, agg_fun = mean)
#' }
#' @export
plotCCIsummary <- function(x, ...) UseMethod("plotCCIsummary")


#' @describeIn plotCCIsummary Method for a \code{blisa} object. Stops with an
#'   informative error if \code{CCI_scores} is \code{NULL}.
#'
#' @param top_lr Integer or \code{NULL}. Number of top-ranked LR pairs (by
#'   \code{sig_numbers}) to include before aggregating. LR pairs in
#'   \code{CCI_scores} are already ordered by rank, so this takes the first
#'   \code{top_lr} columns. \code{NULL} (default) uses all pairs.
#' @param sender Character vector or \code{NULL}. If provided, only rows where
#'   \code{Sender} is in this vector are kept (AND logic with \code{receiver}).
#'   Default \code{NULL} (all senders).
#' @param receiver Character vector or \code{NULL}. If provided, only rows
#'   where \code{Receiver} is in this vector are kept (AND logic with
#'   \code{sender}). Default \code{NULL} (all receivers).
#' @param agg_fun Function used to aggregate scores across LR pairs for each
#'   Sender \eqn{\rightarrow} Receiver combination. Receives a numeric vector
#'   with \code{NA}s already removed. Default \code{sum}.
#' @param main Character or \code{NULL}. Title drawn above the heatmap. When
#'   supplied, the heatmap is drawn with this overall title (via
#'   \code{ComplexHeatmap::draw}); the \code{Heatmap} object is returned
#'   invisibly. Default \code{NULL} (no title; object returned for the caller
#'   to print).
#'
#' @export
plotCCIsummary.blisa <- function(x, top_lr = NULL, sender = NULL,
                                  receiver = NULL, agg_fun = sum, main = NULL,
                                  ...) {
  if (is.null(x$CCI_scores))
    stop("CCI_scores is NULL. Run runCCI() first to compute CCI scores.")
  plotCCIsummary.default(x$CCI_scores, top_lr = top_lr, sender = sender,
                          receiver = receiver, agg_fun = agg_fun, main = main)
}


#' @describeIn plotCCIsummary Method for a CCI scores data frame (e.g. the
#'   \code{CCI_scores} slot of a \code{blisa} object).
#'
#' @export
plotCCIsummary.default <- function(x, top_lr = NULL, sender = NULL,
                                    receiver = NULL, agg_fun = sum, main = NULL,
                                    ...) {
  CCI_df <- x

  if (!is.null(sender)) {
    CCI_df <- CCI_df[CCI_df$Sender %in% sender, , drop = FALSE]
    if (nrow(CCI_df) == 0)
      stop("No rows remaining after filtering by sender.")
  }

  if (!is.null(receiver)) {
    CCI_df <- CCI_df[CCI_df$Receiver %in% receiver, , drop = FALSE]
    if (nrow(CCI_df) == 0)
      stop("No rows remaining after filtering by receiver.")
  }

  lr_cols <- setdiff(colnames(CCI_df), c("Sender", "Receiver"))
  if (!is.null(top_lr))
    lr_cols <- lr_cols[seq_len(min(top_lr, length(lr_cols)))]

  score_mat <- as.matrix(CCI_df[, lr_cols, drop = FALSE])
  agg_scores <- apply(score_mat, 1, function(v) agg_fun(v[!is.na(v)]))

  senders   <- CCI_df$Sender
  receivers <- CCI_df$Receiver
  all_senders   <- unique(senders)
  all_receivers <- unique(receivers)

  interaction_mat <- matrix(NA,
                            nrow = length(all_receivers),
                            ncol = length(all_senders),
                            dimnames = list(all_receivers, all_senders))

  for (i in seq_along(agg_scores)) {
    interaction_mat[receivers[i], senders[i]] <- agg_scores[i]
  }

  f1 <- viridisLite::viridis(10)
  p <- Heatmap(interaction_mat,
               name = "Interaction\nScore",
               col = f1,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               row_title = "Receiver",
               column_title = "Sender",
               row_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_names_rot = 45,
               na_col = "grey90",
               heatmap_legend_param = list(title_position = "topcenter"))

  if (!is.null(main)) {
    draw(p, column_title = main,
         column_title_gp = gpar(fontsize = 14, fontface = "bold"))
    return(invisible(p))
  }
  return(p)
}
