#' Heatmap of CCI scores for a single ligand-receptor pair
#'
#' Reshapes the CCI data frame into a receiver-by-sender matrix for one LR pair
#' and draws a clustered heatmap.
#'
#' @param CCIres The \code{CCI_scores} slot of a \code{blisa} object (i.e.
#'   \code{res$CCI_scores}). Must contain columns \code{Sender}, \code{Receiver},
#'   and one column per LR pair.
#' @param lr_pair Character. Column name in \code{CCIres} corresponding to the
#'   ligand-receptor pair to visualise (e.g. \code{"CD80_CD274"}).
#'
#' @return A \code{Heatmap} object.
#' @export
CCIheatmapOneLR <- function(CCIres, lr_pair) {
  # CCIres: your pre-calculated dataframe (cell-cell pairs as rows, LR pairs as columns)
  # lr_pair: column name of the LR pair to visualize, e.g. "CD80_CD274"

  # Check column
  if (!lr_pair %in% colnames(CCIres)) {
    stop(paste("LR pair", lr_pair, "not found in dataframe."))
  }

  # Extract interaction scores for the selected LR pair
  interaction_scores <- CCIres[[lr_pair]]

  sender   <- CCIres$Sender
  receiver <- CCIres$Receiver

  # Create matrix (receiver as rows, sender as columns)
  all_receivers <- unique(receiver)
  all_senders <- unique(sender)

  interaction_mat <- matrix(NA,
                            nrow = length(all_receivers),
                            ncol = length(all_senders),
                            dimnames = list(all_receivers, all_senders))

  for (i in seq_along(interaction_scores)) {
    r <- receiver[i]
    s <- sender[i]
    interaction_mat[r, s] <- interaction_scores[i]
  }

  # Define color scale
  f1 = viridisLite::viridis(10)
  # Plot
  p <- Heatmap(interaction_mat,
               name = "Interaction\nScore",
               col = f1,
               cluster_rows = TRUE,
               cluster_columns = TRUE,
               row_title = paste0(lr_pair, " - Receiver"),
               column_title = paste0(lr_pair, " - Sender"),
               row_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_names_gp = gpar(fontsize = 12, fontface = "bold"),
               column_names_rot = 45,
               heatmap_legend_param = list(title_position = "topcenter"))

  return(p)
}
