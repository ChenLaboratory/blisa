#' Spatial map of dominant sender-receiver cell-type pairs at BLISA hotspots
#'
#' For a selected ligand-receptor pair, identifies the dominant interacting
#' cell-type pair at each significant hotspot bin and draws a spatial map of
#' the tissue coloured by those pairs. Receiver cells are those inside hotspot
#' bins; sender cells are drawn from the immediate neighbourhood.
#'
#' @param x A \code{blisa} object as returned by \code{\link{blisa}}.
#' @param spe A cell-level \code{SpatialExperiment} object.
#' @param index Integer. Row index into \code{x$LR_results} selecting the
#'   ligand-receptor pair to visualise. Ignored when both \code{ligand} and
#'   \code{receptor} are supplied. Default \code{1} (top-ranked pair).
#' @param ligand Character. Ligand gene symbol. When both \code{ligand} and
#'   \code{receptor} are provided the matching LR pair is located automatically
#'   and \code{index} is ignored. Must be supplied together with
#'   \code{receptor}.
#' @param receptor Character. Receptor gene symbol. Must be supplied together
#'   with \code{ligand}.
#' @param ct_group Character. Column name in \code{colData(spe)} containing
#'   cell-type labels. Default \code{"cell_type"}.
#' @param top Integer. Maximum number of distinct cell-type pairs to show in
#'   the legend; remaining pairs are grouped as \code{"rare pairs"}. Default
#'   \code{30}.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotHotspots}} for a significance-based spatial map of
#'   hotspot bins.
#' @export
plotCCIspatial <- function(x, spe, index = 1, ligand = NULL, receptor = NULL,
                            ct_group = "cell_type", top = 30) {

  LR_results <- x$LR_results
  bins       <- x$bins
  sw         <- x$spatial_weights

  index   <- .resolve_lr_index(LR_results, index, ligand, receptor)
  lr_pair <- rownames(LR_results)[index]

  genes <- unname(unlist(LR_results[index, c("ligand.symbol", "receptor.symbol")]))
  sigHH <- LR_results$sig_index[[index]]
  mode  <- LR_results$ccc_mode[index]

  if (length(sigHH) == 0)
    stop("LR pair '", lr_pair, "' has no significant hotspots.")

  nb_list <- if (mode == "nearby") sw$queen_nb_full else sw$dist_nb_full

  # Map each cell to its bin row position
  cell_to_hex <- get_cell_hex_mapping(spe, bins)

  cell_data <- data.frame(
    hex_id        = as.integer(cell_to_hex),
    ct            = as.character(SummarizedExperiment::colData(spe)[[ct_group]]),
    ligand_expr   = as.numeric(SummarizedExperiment::assay(spe, "counts")[genes[1], ]),
    receptor_expr = as.numeric(SummarizedExperiment::assay(spe, "counts")[genes[2], ])
  )

  # Receiver: sum receptor expression per (hotspot bin, cell type)
  rcpt_data    <- cell_data[cell_data$hex_id %in% sigHH, ]
  rcpt_summary <- aggregate(receptor_expr ~ hex_id + ct,
                            data = rcpt_data, FUN = sum)
  colnames(rcpt_summary)[colnames(rcpt_summary) == "ct"]            <- "ct_r"
  colnames(rcpt_summary)[colnames(rcpt_summary) == "receptor_expr"] <- "r_sum"

  # Sender: sum ligand expression per (hotspot + neighbour bin, cell type)
  sigHH_ng    <- unique(c(sigHH, unlist(nb_list[sigHH])))
  lig_data    <- cell_data[cell_data$hex_id %in% sigHH_ng, ]
  lig_summary <- aggregate(ligand_expr ~ hex_id + ct,
                           data = lig_data, FUN = sum)
  colnames(lig_summary)[colnames(lig_summary) == "ct"]           <- "ct_l"
  colnames(lig_summary)[colnames(lig_summary) == "hex_id"]       <- "hh_hex"
  colnames(lig_summary)[colnames(lig_summary) == "ligand_expr"]  <- "l_sum"

  # Cartesian merge: all (receiver ct x sender ct) combinations per hotspot bin
  merged           <- merge(rcpt_summary, lig_summary,
                            by.x = "hex_id", by.y = "hh_hex")
  merged$product   <- 0.5 * log2(merged$r_sum * merged$l_sum + 1)
  merged$cell_pair <- paste(merged$ct_l, merged$ct_r, sep = " -> ")

  # Keep the top-scoring sender-receiver pair per hotspot bin
  best_idx  <- tapply(seq_len(nrow(merged)), merged$hex_id,
                      function(i) i[which.max(merged$product[i])])
  top_pairs <- merged[unlist(best_idx), ]

  # Legend: show top N pairs by frequency; group remainder as "rare pairs"
  tbl <- sort(table(top_pairs$cell_pair), decreasing = TRUE)
  if (length(tbl) <= top) {
    shown_pairs  <- names(tbl)
    legend_title <- "All pairs"
  } else {
    shown_pairs  <- names(tbl[seq_len(top)])
    legend_title <- paste0("Top ", top, " pairs")
  }
  top_pairs$cell_pair_plot <- ifelse(top_pairs$cell_pair %in% shown_pairs,
                                     top_pairs$cell_pair, "rare pairs")

  # Label every bin: Empty / Non-Significant / dominant interacting pair
  bins$cell_pair_plot                    <- ifelse(bins$n_cells > 0,
                                                    "Non-Significant", "Empty")
  bins$cell_pair_plot[top_pairs$hex_id] <- top_pairs$cell_pair_plot
  legend_levels       <- c(shown_pairs, "rare pairs", "Non-Significant", "Empty")
  bins$cell_pair_plot <- factor(bins$cell_pair_plot, levels = legend_levels)

  fill_values <- c(
    setNames(cols[seq_along(shown_pairs)], shown_pairs),
    "rare pairs"      = "#818589",
    "Non-Significant" = "#D3D3D3",
    "Empty"           = "#F0F0F0"
  )

  p <- ggplot(bins) +
    geom_sf(aes(fill = cell_pair_plot), color = NA) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    guides(fill = guide_legend(title = legend_title)) +
    labs(
      title    = paste0(lr_pair, ": Interacting Hotspots"),
      subtitle = "Grey: light = empty bins, medium = non-significant bins"
    ) +
    theme_void()

  return(p)
}
