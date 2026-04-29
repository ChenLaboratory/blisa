get_cell_hex_mapping <- function(spe, hex_sf) {
  # Convert spe coordinates to sf
  coords <- as.data.frame(SpatialExperiment::spatialCoords(spe))
  cell_sf <- sf::st_as_sf(coords, coords = c("x_centroid", "y_centroid"), crs = sf::st_crs(hex_sf))

  # Spatial join to find which hex each cell falls into
  # join = st_intersects is the standard for points in polygons
  mapped <- sf::st_join(cell_sf, hex_sf, join = sf::st_intersects)

  # Return a named vector: names = cell_id, values = hex_id
  res <- mapped$hex_id
  names(res) <- rownames(coords)
  return(res)
}

LR_source_df_updated <- function(spe, LRI_sum, cell_to_hex, index, queen_nb, dist_nb, ct_group = "cell_type") {

  # 1. Get Genes and Mode
  l_gene <- LRI_sum$ligand.symbol[index]
  r_gene <- LRI_sum$receptor.symbol[index]
  mode   <- LRI_sum$ccc_mode[index]
  sigHH  <- LRI_sum$sig_index[[index]] # Hotspot indices

  # 2. Identify Neighbors based on Mode
  # Use the same neighbor lists (nb) used in runBLISA
  nb_list <- if (mode == "nearby") queen_nb else dist_nb

  # Get all neighbors of the hotspots
  # unique(unlist(nb_list[sigHH])) finds all hexagons within the defined range
  sigHH_ng <- unique(c(sigHH, unlist(nb_list[sigHH])))
  # Note: we include sigHH in the sender list because cells can talk to themselves/neighbors

  # 3. Aggregate counts by Cell Type
  # Receptor/Receiver: Only cells inside the Hotspots
  cells_in_HH <- names(cell_to_hex)[cell_to_hex %in% sigHH]

  # Ligand/Sender: Cells in the Hotspots AND their defined neighbors
  cells_in_HH_ng <- names(cell_to_hex)[cell_to_hex %in% sigHH_ng]

  # Helper to sum counts
  get_sums <- function(target_cells, gene) {
    if(length(target_cells) == 0) return(numeric(0))
    # Using aggregate for speed on Xenium cell-level data
    counts_vec <- as.numeric(counts(spe)[gene, target_cells])
    cts_vec <- SummarizedExperiment::colData(spe)[target_cells, ct_group]
    tapply(counts_vec, cts_vec, sum)
  }

  rec_sums <- get_sums(cells_in_HH, r_gene)
  send_sums <- get_sums(cells_in_HH_ng, l_gene)

  # 4. Final Formatting
  all_cts <- unique(SummarizedExperiment::colData(spe)[[ct_group]])
  res_df <- data.frame(
    receiver_sum = as.numeric(rec_sums[all_cts]),
    sender_sum = as.numeric(send_sums[all_cts]),
    row.names = all_cts
  )
  res_df[is.na(res_df)] <- 0

  return(res_df)
}

#' Score cell-cell interactions from BLISA hotspots
#'
#' For each significant ligand-receptor pair, aggregates ligand expression in
#' sender cells (hotspot neighbourhood) and receptor expression in receiver cells
#' (hotspot bins) by cell type, then computes a geometric-mean interaction score.
#' Returns a wide data frame: rows are \code{"Sender->Receiver"} cell-type pairs,
#' columns are LR pairs.
#'
#' @param spe A cell-level \code{SpatialExperiment} object.
#' @param BLISA_output Result list returned by \code{runBLISA.spe} or
#'   \code{runBLISA.spe.isolates.removed}. Must contain \code{LR_out} and
#'   \code{hex_sf}.
#' @param ct_group Character. Column name in \code{colData(spe)} containing
#'   cell-type labels. Default \code{"cell_type"}.
#' @param hex_size Numeric. Hex bin spacing used to compute queen neighbours for
#'   nearby-mode interactions. Default 50.
#' @param dmax Numeric. Maximum distance used to compute distance neighbours for
#'   diffuse-mode interactions. Default 250.
#'
#' @return A data frame with \code{"Sender->Receiver"} row names and one column
#'   per LR pair containing the interaction score.
#' @export
runCCI <- function(spe, BLISA_output, ct_group = "cell_type", hex_size = 50, dmax = 250) {

  # Extract components from your BLISA result list
  LRI_sum <- BLISA_output$LR_out
  hex_sf  <- BLISA_output$hex_sf

  # 1. Pre-calculate Neighbors (Exact same logic as runBLISA)
  centroids <- sf::st_centroid(hex_sf)
  coords <- sf::st_coordinates(centroids)

  dist_nb <- spdep::dnearneigh(coords, 0, dmax)
  queen_nb <- spdep::dnearneigh(centroids, 0, 1.2 * hex_size)

  # 2. Map cells to hexagons
  cell_to_hex <- get_cell_hex_mapping(spe, hex_sf)

  # 3. Calculate CCI for all pairs
  interaction_list <- lapply(seq_len(nrow(LRI_sum)), function(idx) {
    if (LRI_sum$sig_numbers[idx] == 0) return(NULL) # Skip pairs with no hotspots

    LR_source_sum <- LR_source_df_updated(spe, LRI_sum, cell_to_hex, idx,
                                          queen_nb, dist_nb, ct_group)

    # Calculate geometric mean of log-transformed sums
    receiver <- log10(LR_source_sum$receiver_sum + 1)
    sender   <- log10(LR_source_sum$sender_sum + 1)

    score_mat <- sqrt(outer(receiver, sender, FUN = "*"))

    rownames(score_mat) <- rownames(LR_source_sum)  # Receiver
    colnames(score_mat) <- rownames(LR_source_sum)  # Sender

    df <- as.data.frame(as.table(score_mat))
    colnames(df) <- c("Receiver", "Sender", "Score")
    df$LR_pair <- rownames(LRI_sum)[idx]
    df$CellPair <- paste(df$Sender, df$Receiver, sep = "->")
    return(df)
  })

  # Combine and Reshape
  all_df <- do.call(rbind, interaction_list)
  interaction_wide <- reshape2::dcast(all_df, CellPair ~ LR_pair, value.var = "Score", fill = 0)
  rownames(interaction_wide) <- interaction_wide$CellPair

  return(interaction_wide[, -1, drop = FALSE])
}
