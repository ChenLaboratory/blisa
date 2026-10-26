#' Spatial map of dominant sender-receiver cell-type pairs at BLISA hotspots
#'
#' For a selected ligand-receptor pair, identifies the dominant interacting
#' cell-type pair at each significant hotspot bin and draws a spatial map of
#' the tissue coloured by those pairs. Receiver cells are those inside hotspot
#' bins; sender cells are drawn from the neighbouring bins defined by the
#' LR pair's signalling mode.
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
#' @param weighted Logical. If \code{TRUE} (default), the ligand contribution
#'   of each sender cell type at hotspot \eqn{h} is
#'   \eqn{\sum_{j \in N(h)} w_{hj} L_j}, using the same neighbours and
#'   row-standardised spatial weights as the bivariate Moran's I test (the
#'   hotspot bin itself is excluded). If \code{FALSE}, ligand counts are summed
#'   without weights over the hotspot bin and its neighbours.
#'
#' @return A \code{ggplot} object.
#' @seealso \code{\link{plotHotspots}} for a significance-based spatial map of
#'   hotspot bins.
#' @examples
#' \dontrun{
#' # Continuing from the blisa() example:
#' # result <- blisa(spe, bin_size = 50, group = "cell_type")
#' binned <- hexBinCells(
#'   as.data.frame(SpatialExperiment::spatialCoords(spe)),
#'   SummarizedExperiment::assay(spe, "counts"),
#'   bin_size = 50, group = spe$cell_type
#' )
#' plotCCIspatial(result, binned$counts_by_group, index = 1)
#' }
#' @export
plotCCIspatial <- function(x, counts_by_group, index = 1, ligand = NULL,
                            receptor = NULL, top_pairs = 30, weighted = TRUE) {

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

  if (mode == "nearby") {
    nb_list  <- sw$queen_nb_full
    wt       <- sw$queen_wt
    keep_idx <- sw$keep_idx_queen
  } else {
    nb_list  <- sw$dist_nb_full
    wt       <- sw$dist_wt
    keep_idx <- sw$keep_idx_dist
  }

  # Collapse each cell type's count matrix to a single per-bin expression vector
  # ONCE (get_min_expr takes the per-bin minimum across multi-subunit complexes).
  # Doing this outside the hotspot loop avoids recomputing the full-width vector
  # for every hotspot bin.
  r_expr <- lapply(counts_by_group, function(m) get_min_expr(gene_r, m))
  l_expr <- lapply(counts_by_group, function(m) get_min_expr(gene_l, m))

  # For each hotspot bin, score all sender-receiver cell-type combinations
  # and retain the dominant (highest-scoring) pair
  dominant_pairs <- do.call(rbind, lapply(sigHH, function(h) {
    if (weighted) {
      k    <- match(h, keep_idx)
      nb_h <- keep_idx[wt$neighbours[[k]]]
      w_h  <- wt$weights[[k]]
    } else {
      nb_h <- unique(c(h, nb_list[[h]]))
      w_h  <- 1
    }

    # Receptor expression in bin h per cell type (receivers)
    r_scores <- vapply(ct_names, function(ct) sum(r_expr[[ct]][h]),    numeric(1))
    # Ligand expression in neighbourhood of h per cell type (senders)
    l_scores <- vapply(ct_names, function(ct) sum(w_h * l_expr[[ct]][nb_h]), numeric(1))

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

  # Label every bin: Empty / Non-Significant / dominant interacting pair.
  # Bins not included in this LR pair's LISA test (empty, isolated, low-cell,
  # or low total counts) are treated the same as empty bins.
  tested <- !is.na(LR_results$all_quadrant[[index]])
  bins$cell_pair_plot                        <- ifelse(tested,
                                                       "Non-Significant", "Empty")
  bins$cell_pair_plot[dominant_pairs$hex_id] <- dominant_pairs$cell_pair_plot
  legend_levels       <- c(shown_pairs, "rare pairs", "Non-Significant", "Empty")
  bins$cell_pair_plot <- factor(bins$cell_pair_plot, levels = legend_levels)

  fill_values <- c(
    setNames(cols[seq_along(shown_pairs)], shown_pairs),
    "rare pairs"      = "#404040",
    "Non-Significant" = "#D3D3D3",
    "Empty"           = "#FFFFFF"
  )

  p <- ggplot(bins) +
    geom_sf(aes(fill = cell_pair_plot), color = NA) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    guides(fill = guide_legend(title = legend_title)) +
    labs(title = paste0(lr_pair, ": Interacting Hotspots")) +
    theme_void()

  return(p)
}
