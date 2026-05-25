#' Spatial map of dominant sender-receiver cell-type pairs at BLISA hotspots
#'
#' Plots the dominant interacting cell-type pair at each hotspot bin for a
#' selected ligand-receptor interaction. Receiver cells are those inside hotspot
#' bins; sender cells are in the immediate neighbourhood.
#'
#' @param spe A cell-level \code{SpatialExperiment} object.
#' @param BLISA_output An object of class \code{blisa} as returned by
#'   \code{\link{blisa}}.
#' @param index Integer. Row index into \code{BLISA_output$LR_results} selecting
#'   the ligand-receptor pair to visualise.
#' @param ct_group Character. Column name in \code{colData(spe)} containing
#'   cell-type labels. Default \code{"cell_type"}.
#' @param top Integer. Maximum number of distinct cell-type pairs to show in
#'   the legend; remaining pairs are grouped as \code{"rare pairs"}. Default 30.
#'
#' @return A \code{ggplot} object.
#' @export
CCIspatial <- function(
    spe,
    BLISA_output,
    index,
    ct_group = "cell_type",
    top = 30
) {
  # 1. Setup Data
  LRI_sum         <- BLISA_output$LR_results
  bins          <- BLISA_output$bins
  sw <- BLISA_output$spatial_weights

  interaction <- unname(unlist(LRI_sum[index, c("ligand.symbol", "receptor.symbol")]))
  sigHH <- LRI_sum$sig_index[[index]]
  mode  <- LRI_sum$ccc_mode[index]

  # 2. Neighbour lookup (reuse pre-computed lists from blisa)
  nb_list <- if (mode == "nearby") sw$queen_nb_full else sw$dist_nb_full

  # 3. Map cells to bin row positions
  cell_to_hex <- get_cell_hex_mapping(spe, bins)

  cell_data <- data.table::data.table(
    hex_id = as.integer(cell_to_hex),
    ct = as.character(SummarizedExperiment::colData(spe)[[ct_group]]),
    ligand_expr = as.numeric(SummarizedExperiment::counts(spe)[interaction[1], ]),
    receptor_expr = as.numeric(SummarizedExperiment::counts(spe)[interaction[2], ])
  )

  # 4. Interaction Scoring (Receiver in HH bins, Sender in neighbours)
  rcpt_summary <- cell_data[hex_id %in% sigHH, .(r_sum = sum(receptor_expr)), by = .(hex_id, ct_r = ct)]

  sigHH_ng <- unique(c(sigHH, unlist(nb_list[sigHH])))
  lig_summary <- cell_data[hex_id %in% sigHH_ng, .(l_sum = sum(ligand_expr)), by = .(hh_hex = hex_id, ct_l = ct)]

  # Merge receiver and sender summaries and find top pair per hotspot
  merged_scores <- merge(rcpt_summary, lig_summary, by.x = "hex_id", by.y = "hh_hex", allow.cartesian = TRUE)
  merged_scores[, product := 0.5 * log2(r_sum * l_sum + 1)]
  merged_scores[, cell_pair := paste(ct_l, ct_r, sep = " -> ")]

  top_pairs <- merged_scores[merged_scores[, .I[which.max(product)], by = hex_id]$V1]

  # 5. Legend and Filtering Logic
  tbl <- sort(table(top_pairs$cell_pair), decreasing = TRUE)

  if (length(tbl) <= top) {
    filtered_pairs <- names(tbl)
    legend_title <- "All pairs (ordered)"
  } else {
    filtered_pairs <- names(tbl[1:top])
    legend_title <- paste0("Top ", top, " pairs (ordered)")
  }

  # label pairs outside the top-N as "rare pairs" for legend clarity
  top_pairs[, cell_pair_plot := ifelse(cell_pair %in% filtered_pairs, cell_pair, "rare pairs")]

  # 6. Categorize ALL bins
  bins$cell_pair_plot <- ifelse(bins$n_cells > 0, "Non-Significant", "Empty")

  # Overwrite hotspot bins with their dominant interacting pair
  # top_pairs$hex_id are row positions in bins (from get_cell_hex_mapping)
  bins$cell_pair_plot[top_pairs$hex_id] <- top_pairs$cell_pair_plot

  # Order levels for the legend
  legend_levels <- c(filtered_pairs, "rare pairs", "Non-Significant", "Empty")
  bins$cell_pair_plot <- factor(bins$cell_pair_plot, levels = legend_levels)

  # 7. Visualization
  p <- ggplot(bins) +
    geom_sf(aes(fill = cell_pair_plot), color = NA) +
    scale_fill_manual(
      values = c(
        "Empty" = "#F0F0F0",
        "Non-Significant" = "#D3D3D3",
        "rare pairs" = "#818589",
        setNames(cols[seq_along(filtered_pairs)], filtered_pairs)
      )
    ) +
    theme_void() +
    guides(fill = guide_legend(title = legend_title)) +
    labs(
      title = paste0(rownames(LRI_sum)[index], ": Interacting Hotspots"),
      subtitle = "Grey scale: Light (Empty) | Medium (Non-sig cells)"
    )

  return(p)
}
