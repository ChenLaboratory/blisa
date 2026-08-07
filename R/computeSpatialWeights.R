#' Compute Spatial Weights for BLISA
#'
#' Builds queen (nearby) and distance-decay (diffuse) spatial weight matrices
#' from a bin-level \code{sf} object, excluding bins that have no neighbours.
#'
#' @param bins An \code{sf} object of spatial bins.
#' @param bin_size Numeric. Bin spacing used to define queen adjacency
#'   (\code{1.2 * bin_size} radius).
#' @param dmax Numeric. Maximum distance for diffuse-mode neighbours.
#'
#' @return A list with:
#' \describe{
#'   \item{queen_wt}{Spatial weights list for nearby (queen) mode.}
#'   \item{dist_wt}{Spatial weights list for diffuse (distance-decay) mode.}
#'   \item{keep_idx_queen}{Integer indices of bins used in queen-mode Moran.}
#'   \item{keep_idx_dist}{Integer indices of bins used in diffuse-mode Moran.}
#'   \item{isolate_idx_queen}{Integer indices of original queen-mode isolates.}
#'   \item{isolate_idx_dist}{Integer indices of original diffuse-mode isolates.}
#'   \item{queen_nb_full}{Full (unsubset) neighbour list for nearby mode, indexed over all bins.}
#'   \item{dist_nb_full}{Full (unsubset) neighbour list for diffuse mode, indexed over all bins.}
#' }
#' @examples
#' \dontrun{
#' set.seed(42)
#' pts  <- sf::st_as_sf(
#'   data.frame(x = runif(300, 0, 1000), y = runif(300, 0, 1000)),
#'   coords = c("x", "y"), crs = NA
#' )
#' bins <- sf::st_sf(
#'   geometry = sf::st_make_grid(pts, cellsize = 100,
#'                               what = "polygons", square = FALSE)
#' )
#' sw <- computeSpatialWeights(bins, bin_size = 100, dmax = 300)
#' names(sw)
#' }
#' @export
computeSpatialWeights <- function(bins,
                                  bin_size = 50,
                                  dmax     = 250) {
  centroids <- sf::st_centroid(bins)
  coords    <- sf::st_coordinates(centroids)
  n_bins    <- nrow(bins)

  ## ---------------------------
  ## Queen spatial weights  (for "nearby" mode)
  ## ---------------------------
  queen_nb_full     <- spdep::dnearneigh(coords, 0, 1.2 * bin_size)
  isolate_idx_queen <- which(spdep::card(queen_nb_full) == 0)
  message(length(isolate_idx_queen), " isolated bins with no nearby neighbours: ",
          paste(isolate_idx_queen, collapse = ","))

  keep_idx_queen <- setdiff(seq_len(n_bins), isolate_idx_queen)
  queen_nb       <- spdep::subset.nb(queen_nb_full,
                                     subset = seq_len(n_bins) %in% keep_idx_queen)

  queen_wt <- spdep::nb2listwdist(queen_nb, centroids[keep_idx_queen, ],
                                  type = "idw", style = "W", zero.policy = TRUE)

  ## ---------------------------
  ## Distance spatial weights  (for "diffuse" mode)
  ## ---------------------------
  dist_nb_full     <- spdep::dnearneigh(coords, 0, dmax)
  isolate_idx_dist <- which(spdep::card(dist_nb_full) == 0)
  message(length(isolate_idx_dist), " isolated bins with no neighbours within ",
          dmax, " um: ", paste(isolate_idx_dist, collapse = ","))

  keep_idx_dist <- setdiff(seq_len(n_bins), isolate_idx_dist)
  dist_nb       <- spdep::subset.nb(dist_nb_full,
                                    subset = seq_len(n_bins) %in% keep_idx_dist)

  weight_at_dmax <- 0.01
  dist_wt <- spdep::nb2listwdist(dist_nb, centroids[keep_idx_dist, ],
                                 type = "exp", style = "W", zero.policy = TRUE,
                                 alpha = -log(weight_at_dmax) / dmax)

  list(
    queen_wt          = queen_wt,
    dist_wt           = dist_wt,
    keep_idx_queen    = keep_idx_queen,
    keep_idx_dist     = keep_idx_dist,
    isolate_idx_queen = isolate_idx_queen,
    isolate_idx_dist  = isolate_idx_dist,
    queen_nb_full     = queen_nb_full,
    dist_nb_full      = dist_nb_full
  )
}
