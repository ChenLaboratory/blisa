get_cell_hex_mapping <- function(spe, bins) {
  # Convert spe coordinates to sf
  coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
  cell_sf <- sf::st_as_sf(coords, coords = c("x_centroid", "y_centroid"), crs = sf::st_crs(bins))

  # Spatial join to find which bin each cell falls into
  mapped <- sf::st_join(cell_sf, bins, join = sf::st_intersects)

  # Return row positions in bins (not bin_id values) so they align with
  # the sigHH indices stored in BLISA_output$LR_results$sig_index
  res <- match(mapped$bin_id, bins$bin_id)
  names(res) <- rownames(coords)
  return(res)
}


#' Score cell-cell interactions from BLISA hotspots
#'
#' For each significant ligand-receptor pair, aggregates ligand expression in
#' sender bins (hotspot neighbourhood) and receptor expression in receiver bins
#' (hotspot bins) by group (e.g. cell type), then computes an interaction score
#' \code{0.5 * log2(receiver * sender + 1)} for every sender-receiver group
#' combination. Returns a wide data frame: rows are \code{"Sender->Receiver"}
#' group pairs, columns are LR pairs.
#'
#' @param BLISA_output An object of class \code{blisa} as returned by
#'   \code{\link{blisa}}.
#' @param counts_by_group Named list of gene-by-bin sparse count matrices, one
#'   per group level (e.g. cell type). Typically the \code{counts_by_group}
#'   element of the list returned by \code{\link{hexBinCells}} when \code{group}
#'   is supplied. Names must match the group levels.
#'
#' @return A data frame with \code{"Sender->Receiver"} row names and one column
#'   per LR pair (only pairs with at least one hotspot) containing the
#'   interaction score.
#' @export
runCCI <- function(BLISA_output, counts_by_group) {

  LRI_sum         <- BLISA_output$LR_results
  sw            <- BLISA_output$spatial_weights
  queen_nb_full <- sw$queen_nb_full
  dist_nb_full  <- sw$dist_nb_full
  ct_names      <- names(counts_by_group)

  scores_list <- lapply(seq_len(nrow(LRI_sum)), function(idx) {
    if (LRI_sum$sig_numbers[idx] == 0) return(NULL)

    l_gene <- LRI_sum$ligand.symbol[idx]
    r_gene <- LRI_sum$receptor.symbol[idx]
    mode   <- LRI_sum$ccc_mode[idx]
    sigHH  <- LRI_sum$sig_index[[idx]]

    nb_list  <- if (mode == "nearby") queen_nb_full else dist_nb_full
    sigHH_ng <- unique(c(sigHH, unlist(nb_list[sigHH])))

    # Sum receptor counts in hotspot bins per group (receivers)
    receiver_sums <- sapply(ct_names, function(ct) {
      sum(counts_by_group[[ct]][r_gene, sigHH])
    })

    # Sum ligand counts in hotspot + neighbour bins per group (senders)
    sender_sums <- sapply(ct_names, function(ct) {
      sum(counts_by_group[[ct]][l_gene, sigHH_ng])
    })

    score_mat  <- 0.5 * log2(outer(receiver_sums, sender_sums, FUN = "*") + 1)
    pair_names <- as.vector(outer(ct_names, ct_names, function(r, s) paste(s, r, sep = "->")))
    scores     <- as.vector(score_mat)
    names(scores) <- pair_names
    scores
  })
  names(scores_list) <- rownames(LRI_sum)

  scores_list <- Filter(Negate(is.null), scores_list)
  if (length(scores_list) == 0) return(data.frame())

  as.data.frame(do.call(cbind, scores_list))
}
