#' Sender-by-receiver heatmap of CCI scores for one ligand-receptor pair
#'
#' Generic function. Reshapes the CCI data frame into a receiver-by-sender
#' cell-type matrix for one selected LR pair and draws a clustered heatmap.
#'
#' @param x A \code{blisa} object or a CCI scores data frame (the
#'   \code{CCI_scores} slot of a \code{blisa} object).
#' @param ... Additional arguments passed to the relevant method.
#'
#' @return A \code{Heatmap} object.
#' @seealso \code{\link{plotCCI}} for an overview heatmap across all LR pairs;
#'   \code{\link{plotCCIsummary}} for an aggregated sender-by-receiver heatmap.
#' @export
plotCCIpair <- function(x, ...) UseMethod("plotCCIpair")


#' @describeIn plotCCIpair Method for a \code{blisa} object. The LR pair is
#'   selected by \code{index} (default 1, the top-ranked pair) unless both
#'   \code{ligand} and \code{receptor} are supplied, in which case the matching
#'   row is located automatically and \code{index} is ignored. Stops with an
#'   informative error if \code{CCI_scores} is \code{NULL} or the selected LR
#'   pair has no significant hotspots.
#'
#' @param index Integer. Row index into \code{LR_results} selecting the
#'   ligand-receptor pair to visualise. Ignored when both \code{ligand} and
#'   \code{receptor} are supplied. Default \code{1} (top-ranked pair).
#' @param ligand Character. Ligand gene symbol. When both \code{ligand} and
#'   \code{receptor} are provided the matching LR pair is located automatically
#'   and \code{index} is ignored. Must be supplied together with
#'   \code{receptor}.
#' @param receptor Character. Receptor gene symbol. Must be supplied together
#'   with \code{ligand}.
#'
#' @export
plotCCIpair.blisa <- function(x, index = 1, ligand = NULL, receptor = NULL,
                               ...) {
  if (is.null(x$CCI_scores))
    stop("CCI_scores is NULL. Run runCCI() first to compute CCI scores.")

  LR_results <- x$LR_results

  # Resolve which LR pair to plot (reuses helper from plotHotspots.R)
  index  <- .resolve_lr_index(LR_results, index, ligand, receptor)
  lr_pair <- rownames(LR_results)[index]

  if (LR_results$sig_numbers[index] == 0)
    stop("LR pair '", lr_pair, "' has no significant hotspots and therefore ",
         "no CCI scores.")

  if (!lr_pair %in% colnames(x$CCI_scores))
    stop("LR pair '", lr_pair, "' not found in CCI_scores.")

  plotCCIpair.default(x$CCI_scores, lr_pair)
}


#' @describeIn plotCCIpair Method for a CCI scores data frame (e.g. the
#'   \code{CCI_scores} slot of a \code{blisa} object). The LR pair is selected
#'   by column name via \code{lr_pair}.
#'
#' @param lr_pair Character. Column name in the CCI scores data frame
#'   corresponding to the ligand-receptor pair to visualise
#'   (e.g. \code{"CXCL12_CXCR4"}).
#'
#' @export
plotCCIpair.default <- function(x, lr_pair, ...) {
  CCI_df <- x

  if (!lr_pair %in% colnames(CCI_df))
    stop("LR pair '", lr_pair, "' not found in CCI_df.")

  interaction_scores <- CCI_df[[lr_pair]]
  sender   <- CCI_df$Sender
  receiver <- CCI_df$Receiver

  # Build receiver-by-sender matrix
  all_receivers <- unique(receiver)
  all_senders   <- unique(sender)

  interaction_mat <- matrix(NA,
                            nrow = length(all_receivers),
                            ncol = length(all_senders),
                            dimnames = list(all_receivers, all_senders))

  for (i in seq_along(interaction_scores)) {
    interaction_mat[receiver[i], sender[i]] <- interaction_scores[i]
  }

  f1 <- viridisLite::viridis(10)
  p <- Heatmap(interaction_mat,
               name = "Interaction\nScore",
               col = f1,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               row_title = paste0(lr_pair, " — Receiver"),
               column_title = paste0(lr_pair, " — Sender"),
               row_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_names_rot = 45,
               heatmap_legend_param = list(title_position = "topcenter"))

  return(p)
}
