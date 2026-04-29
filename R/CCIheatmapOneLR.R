
# Define global color palette
# heatmap_cols <- viridisLite::viridis(256, option = "viridis")
#
# CCIheatmapOneLR <- function(
#     CCIres,
#     LRI
# ) {
#   # 1. Validation: Check if the LRI exists in your results
#   if (!(LRI %in% colnames(CCIres))) {
#     stop(paste("LRI", LRI, "not found in the CCI matrix. Available pairs include:",
#                paste(head(colnames(CCIres)), collapse = ", ")))
#   }
#
#   message("Generating heatmap for interaction: ", LRI)
#
#   # 2. Extract and Transform Data
#   # Convert the single LRI column into a data frame with Sender/Receiver columns
#   df_long <- data.frame(
#     CellPair = rownames(CCIres),
#     Score = CCIres[, LRI]
#   )
#
#   # Split "Sender->Receiver" into two distinct variables
#   df_split <- df_long %>%
#     tidyr::separate(CellPair, into = c("Sender", "Receiver"), sep = "->")
#
#   # Pivot into a wide matrix: Rows = Receiver, Cols = Sender
#   # Use acast for a clean matrix output for ComplexHeatmap
#   interaction_mat <- reshape2::acast(df_split, Receiver ~ Sender, value.var = "Score")
#
#   # 3. Dynamic Color Mapping
#   # This ensures the 256 viridis levels are mapped to the actual range of this LRI
#   val_range <- range(interaction_mat, na.rm = TRUE)
#   if (val_range[1] == val_range[2]) {
#     val_range[2] <- val_range[1] + 1 # Prevent colorRamp2 error on constant data
#   }
#
#   col_fun = colorRamp2(
#     seq(val_range[1], val_range[2], length.out = 256),
#     heatmap_cols
#   )
#
#   # 4. Parsing Gene Names for Titles
#   # Assumes format "Ligand_Receptor"
#   genes <- unlist(strsplit(LRI, "_"))
#   l_gene <- genes[1]
#   r_gene <- genes[2]
#
#   # 5. Build the Heatmap
#   p <- Heatmap(
#     interaction_mat,
#     name = "Interaction\nScore",
#     col = col_fun,
#     cluster_rows = TRUE,
#     cluster_columns = TRUE,
#     row_title = paste0(r_gene, " (Receiver Cell Type)"),
#     column_title = paste0(l_gene, " (Sender Cell Type)"),
#     row_names_gp = gpar(fontsize = 11, fontface = "bold"),
#     column_names_gp = gpar(fontsize = 11, fontface = "bold"),
#     column_names_rot = 45,
#     heatmap_legend_param = list(
#       title_position = "topcenter",
#       legend_width = unit(4, "cm"),
#       direction = "vertical"
#     )
#   )
#
#   return(p)
# }

#' Heatmap of CCI scores for a single ligand-receptor pair
#'
#' Reshapes the CCI data frame into a receiver-by-sender matrix for one LR pair
#' and draws a clustered heatmap.
#'
#' @param CCIres Data frame returned by \code{runCCI}. Rows are
#'   \code{"Sender->Receiver"} cell-type pairs; columns are LR pairs.
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

  # Parse receiver and sender names from rownames (e.g., "B_Cells->T_Cells")
  pairs <- strsplit(rownames(CCIres), "->")
  receiver <- sapply(pairs, `[`, 2)
  sender <- sapply(pairs, `[`, 1)

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
