#' Compute Spatial Weights for BLISA
#'
#' Builds queen (nearby) and distance-decay (diffuse) spatial weight matrices
#' from a bin-level \code{sf} object, excluding isolated bins and optionally
#' excluding low-cell bins. A second-pass isolation check further removes bins
#' that become isolated after the initial subset.
#'
#' @param bin_sf An \code{sf} object of spatial bins.
#' @param hex_size Numeric. Bin spacing used to define queen adjacency
#'   (\code{1.2 * hex_size} radius).
#' @param dmax Numeric. Maximum distance for diffuse-mode neighbours.
#' @param min_cells_per_bin Integer. Minimum cell count for a bin to be
#'   included. Ignored when \code{n_cells_col = NA}.
#' @param n_cells_col Character or \code{NA}. Column name in \code{bin_sf}
#'   holding per-bin cell counts. Set to \code{NA} to skip cell-count
#'   filtering (default).
#'
#' @return A list with:
#' \describe{
#'   \item{queen_wt}{Spatial weights list for nearby (queen) mode.}
#'   \item{dist_wt}{Spatial weights list for diffuse (distance-decay) mode.}
#'   \item{keep_idx_queen}{Integer indices of bins used in queen-mode Moran.}
#'   \item{keep_idx_dist}{Integer indices of bins used in diffuse-mode Moran.}
#'   \item{isolate_idx_queen}{Integer indices of original queen-mode isolates.}
#'   \item{isolate_idx_dist}{Integer indices of original diffuse-mode isolates.}
#'   \item{low_cell_idx}{Integer indices of bins excluded for low cell counts.}
#'   \item{queen_nb_full}{Full (unsubset) neighbour list for nearby mode, indexed over all bins.}
#'   \item{dist_nb_full}{Full (unsubset) neighbour list for diffuse mode, indexed over all bins.}
#' }
#' @export
computeSpatialWeights <- function(bin_sf,
                                  hex_size          = 50,
                                  dmax              = 250,
                                  min_cells_per_bin = 1,
                                  n_cells_col       = NA) {
  centroids <- sf::st_centroid(bin_sf)
  coords    <- sf::st_coordinates(centroids)
  n_bins    <- nrow(bin_sf)

  ## ---------------------------
  ## Filter low-cell bins
  ## ---------------------------
  if (!is.na(n_cells_col)) {
    if (!n_cells_col %in% colnames(bin_sf))
      stop("Column '", n_cells_col, "' not found in bin_sf.")
    low_cell_idx <- which(bin_sf[[n_cells_col]] < min_cells_per_bin)
    message(length(low_cell_idx), " bins removed: < ", min_cells_per_bin,
            " cells (column: '", n_cells_col, "').")
  } else {
    low_cell_idx <- integer(0)
    message("n_cells_col = NA — cell-count filtering skipped.")
  }

  ## ---------------------------
  ## Helper: second-pass isolation check
  ## After subset.nb, bins whose neighbours were all excluded become new isolates.
  ## Detect them and further shrink keep_idx before building weights.
  ## ---------------------------
  resolve_new_isolates <- function(nb, nb_full, keep_idx) {
    new_iso_sub <- which(spdep::card(nb) == 0) # positions in subset space
    if (length(new_iso_sub) == 0) return(list(nb = nb, keep_idx = keep_idx))
    new_iso_full <- keep_idx[new_iso_sub] # back to full-space indices
    message(length(new_iso_full), " bins became isolated after subset — excluded.")
    keep_idx <- setdiff(keep_idx, new_iso_full)
    nb       <- spdep::subset.nb(nb_full, subset = seq_len(n_bins) %in% keep_idx)
    list(nb = nb, keep_idx = keep_idx)
  }

  ## ---------------------------
  ## Queen spatial weights  (for "nearby" mode)
  ## ---------------------------
  queen_nb_full     <- spdep::dnearneigh(coords, 0, 1.2 * hex_size)
  isolate_idx_queen <- which(spdep::card(queen_nb_full) == 0)
  message(length(isolate_idx_queen), " isolated bins with no nearby neighbours: ",
          paste(isolate_idx_queen, collapse = ","))

  keep_idx_queen <- setdiff(seq_len(n_bins), union(isolate_idx_queen, low_cell_idx))
  queen_nb       <- spdep::subset.nb(queen_nb_full,
                                     subset = seq_len(n_bins) %in% keep_idx_queen)
  r              <- resolve_new_isolates(queen_nb, queen_nb_full, keep_idx_queen)
  queen_nb       <- r$nb;  keep_idx_queen <- r$keep_idx

  queen_wt <- spdep::nb2listwdist(queen_nb, centroids[keep_idx_queen, ],
                                  type = "idw", style = "W", zero.policy = TRUE)

  ## ---------------------------
  ## Distance spatial weights  (for "diffuse" mode)
  ## ---------------------------
  dist_nb_full     <- spdep::dnearneigh(coords, 0, dmax)
  isolate_idx_dist <- which(spdep::card(dist_nb_full) == 0)
  message(length(isolate_idx_dist), " isolated bins with no neighbours within ",
          dmax, " um: ", paste(isolate_idx_dist, collapse = ","))

  keep_idx_dist <- setdiff(seq_len(n_bins), union(isolate_idx_dist, low_cell_idx))
  dist_nb       <- spdep::subset.nb(dist_nb_full,
                                    subset = seq_len(n_bins) %in% keep_idx_dist)
  r             <- resolve_new_isolates(dist_nb, dist_nb_full, keep_idx_dist)
  dist_nb       <- r$nb;  keep_idx_dist <- r$keep_idx

  weight_at_dmax <- 0.01
  dist_wt <- spdep::nb2listwdist(dist_nb, bin_sf[keep_idx_dist, ],
                                 type = "exp", style = "W", zero.policy = TRUE,
                                 alpha = -log(weight_at_dmax) / dmax)

  list(
    queen_wt          = queen_wt,
    dist_wt           = dist_wt,
    keep_idx_queen    = keep_idx_queen,
    keep_idx_dist     = keep_idx_dist,
    isolate_idx_queen = isolate_idx_queen,
    isolate_idx_dist  = isolate_idx_dist,
    low_cell_idx      = low_cell_idx,
    queen_nb_full     = queen_nb_full,
    dist_nb_full      = dist_nb_full
  )
}
