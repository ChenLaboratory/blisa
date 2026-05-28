#' Spatial map of dominant sender-receiver cell-type pairs at BLISA hotspots
#'
#' For a selected ligand-receptor pair, identifies the dominant interacting
#' cell-type pair at each significant hotspot bin and draws a spatial map of
#' the tissue coloured by those pairs. Receiver cells are those inside hotspot
#' bins; sender cells are drawn from the immediate neighbourhood.
#'
#' @param x A \code{blisa} object as returned by \code{\link{blisa}}.
#' @param counts_by_group Named list of gene-by-bin count matrices, one per
#'   cell type. Typically the \code{counts_by_group} element returned by
#'   \code{\link{hexBinCells}}. Names must match the cell-type levels.
#' @param index Integer. Row index into \code{x$LR_results} selecting the
#'   ligand-receptor pair to visualise. Ignored when both \code{ligand} and
#'   \code{receptor} are supplied. Default \code{1} (top-ranked pair).
#' @param ligand Character. Ligand gene symbol. When both \code{ligand} and
#'   \code{receptor} are provided the matching LR pair is located automatically
#'   and \code{index} is ignored. Must be supplied together with
#'   \code{receptor}.
#' @param receptor Character. Receptor gene symbol. Must be supplied together
#'   with \code{ligand}.
#' @param top_pairs Integer. Maximum number of distinct cell-type pairs to show
#'   in the legend; remaining pairs are grouped as \code{"rare pairs"}. Default
#'   \code{30}.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotHotspots}} for a significance-based spatial map of
#'   hotspot bins.
#' @export
plotCCIspatial <- function(x, counts_by_group, index = 1, ligand = NULL,
                            receptor = NULL, top_pairs = 30) {

  LR_results <- x$LR_results
  bins       <- x$bins
  sw         <- x$spatial_weights
  ct_names   <- names(counts_by_group)

  index   <- .resolve_lr_index(LR_results, index, ligand, receptor)
  lr_pair <- rownames(LR_results)[index]

  gene_l <- LR_results$ligand.symbol[index]
  gene_r <- LR_results$receptor.symbol[index]
  sigHH  <- LR_results$sig_index[[index]]
  mode   <- LR_results$ccc_mode[index]

  if (length(sigHH) == 0)
    stop("LR pair '", lr_pair, "' has no significant hotspots.")

  nb_list <- if (mode == "nearby") sw$queen_nb_full else sw$dist_nb_full

  # For each hotspot bin, score all sender-receiver cell-type combinations
  # and retain the dominant (highest-scoring) pair
  dominant_pairs <- do.call(rbind, lapply(sigHH, function(h) {
    nb_h <- unique(c(h, nb_list[[h]]))

    # Receptor expression in bin h per cell type (receivers)
    r_scores <- sapply(ct_names, function(ct)
      sum(counts_by_group[[ct]][gene_r, h]))

    # Ligand expression in neighbourhood of h per cell type (senders)
    l_scores <- sapply(ct_names, function(ct)
      sum(counts_by_group[[ct]][gene_l, nb_h]))

    score_mat <- 0.5 * log2(outer(r_scores, l_scores) + 1)
    best      <- which(score_mat == max(score_mat), arr.ind = TRUE)[1L, ]

    data.frame(
      hex_id    = h,
      cell_pair = paste(ct_names[best[2L]], ct_names[best[1L]], sep = " -> "),
      product   = score_mat[best[1L], best[2L]]
    )
  }))

  # Legend: show top N pairs by frequency; group remainder as "rare pairs"
  tbl <- sort(table(dominant_pairs$cell_pair), decreasing = TRUE)
  if (length(tbl) <= top_pairs) {
    shown_pairs  <- names(tbl)
    legend_title <- "All pairs"
  } else {
    shown_pairs  <- names(tbl[seq_len(top_pairs)])
    legend_title <- paste0("Top ", top_pairs, " pairs")
  }
  dominant_pairs$cell_pair_plot <- ifelse(dominant_pairs$cell_pair %in% shown_pairs,
                                     dominant_pairs$cell_pair, "rare pairs")

  # Label every bin: Empty / Non-Significant / dominant interacting pair
  bins$cell_pair_plot                   <- ifelse(bins$n_cells > 0,
                                                   "Non-Significant", "Empty")
  bins$cell_pair_plot[dominant_pairs$hex_id] <- dominant_pairs$cell_pair_plot
  legend_levels       <- c(shown_pairs, "rare pairs", "Non-Significant", "Empty")
  bins$cell_pair_plot <- factor(bins$cell_pair_plot, levels = legend_levels)

  fill_values <- c(
    setNames(cols[seq_along(shown_pairs)], shown_pairs),
    "rare pairs"      = "#404040",
    "Non-Significant" = "#D3D3D3",
    "Empty"           = "#F5F5F5"
  )

  p <- ggplot(bins) +
    geom_sf(aes(fill = cell_pair_plot), color = NA) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    guides(fill = guide_legend(title = legend_title)) +
    labs(title = paste0(lr_pair, ": Interacting Hotspots")) +
    theme_void()

  return(p)
}
