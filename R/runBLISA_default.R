#' Run BLISA analysis on bin level sf
#'
#' @param counts_matrix Gene-by-bin matrix of counts; columns must match `bin_sf` bins.
#' @param bin_sf [sf::sf] object with bin polygons and matching row order to columns of `counts_matrix`.
#' @param LR_df Data frame of ligand–receptor pairs (see Details).
#' @param hex_size Numeric. Hexagon size in microns (same units as coordinates).
#' @param dmax Maximum distance (microns) for distance-based neighbors.
#' @param nsim Number of permutations for significance.
#' @param p_cutoff P-value threshold.
#' @param min_ligand,min_receptor Minimum counts required for ligand/receptor.
#' @param col Column name in `LR_df` for interaction category.
#' @param default_mode Default CCC mode when `col` is missing.
#' @param diffuse_category Character vector of categories treated as “diffuse”.
#'
#' @return A list with components `...` (describe structure).
#' @export
runBLISA.default <- function( # bin-level
  counts_matrix,
  bin_sf,
  LR_df, # use index, remove other info
  hex_size = 50,
  dmax = 250,
  nsim = 999,
  p_cutoff = 0.05,
  min_ligand = 10, min_receptor = 10,
  col = "annotation", default_mode = "diffuse",
  diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
) {

  centroids <- sf::st_centroid(bin_sf)
  coords <- sf::st_coordinates(centroids)

  ## ---------------------------
  ## Queen spatial weights
  ## ---------------------------

  ## Queen neighbors
  queen_nb <- spdep::dnearneigh(coords, 0, 1.2 * hex_size)

  # Isolated bins: no queen neighbors
  is_isolate_queen <- spdep::card(queen_nb) == 0
  isolate_idx_queen <- which(is_isolate_queen)
  message(length(isolate_idx_queen), " isolated bins with no nearby neighbors: ", paste(isolate_idx_queen, collapse = ","))

  # Assign a random neighbor for each isolate
  set.seed(123)  # reproducible
  for (i in isolate_idx_queen) {
    # pick a random bin index excluding itself
    queen_nb[[i]] <- sample(setdiff(seq_along(queen_nb), i), 1)
  }

  # Queen weights (for nearby)
  queen_wt <- spdep::nb2listwdist(queen_nb, centroids, type="idw", style="W", zero.policy = TRUE)

  ## ---------------------------
  ## Distance spatial weights
  ## ---------------------------

  ## Distance neighbors
  dist_nb <- spdep::dnearneigh(coords, 0, dmax)

  # Isolated bins: no queen neighbors
  is_isolate_dist <- spdep::card(dist_nb) == 0
  isolate_idx_dist <- which(is_isolate_dist)
  message(length(isolate_idx_dist), " isolated bins with no neighbors within 250 um: ", paste(isolate_idx_dist, collapse = ","))

  # Assign a random neighbor for each isolate
  set.seed(123)  # reproducible
  for (i in isolate_idx_dist) {
    # pick a random bin index excluding itself
    dist_nb[[i]] <- sample(setdiff(seq_along(dist_nb), i), 1)
  }

  # exponential distance weight
  weight_at_dmax <- 0.01 # exp dist weight=0.01 for dist=dmax

  dist_wt <- spdep::nb2listwdist(
    dist_nb,
    bin_sf,
    type = "exp",
    style = "W",
    alpha = -log(weight_at_dmax) / dmax
  )

  ## ---------------------------
  ## Filter LR pairs
  ## ---------------------------
  # only use LR pairs with at least n counts in at least one bin
  LR_df_filtered <- filterLRpairs(counts = counts_matrix,
                                  min_ligand, min_receptor,
                                  LR_df)

  ## ---------------------------
  ## Local bivariate Moran for all LR pairs
  ## ---------------------------
  LR_out <- LR_df_add_mode(LR_df_filtered, col, default_mode, diffuse_category)

  LR_out$sig_numbers <- integer(nrow(LR_out))
  LR_out$sig_index   <- vector("list", nrow(LR_out))
  LR_out$sig_pval    <- vector("list", nrow(LR_out))
  LR_out$all_pval    <- vector("list", nrow(LR_out))
  LR_out$all_lisa    <- vector("list", nrow(LR_out))

  for (i in seq_len(nrow(LR_out))) {
    message(rownames(LR_out)[i])

    ligand   <- LR_out$ligand.symbol[i]
    receptor <- LR_out$receptor.symbol[i]

    mode <- LR_out$ccc_mode[i]
    wt <- if (mode == "nearby") queen_wt else dist_wt
    isolate_idx <- if (mode == "nearby") isolate_idx_queen else isolate_idx_dist
    message("ccc mode is ", mode)

    # bivariate vectors per hex (min for multi-unit)
    x <- get_min_expr(receptor, counts_matrix)  # receptor
    y <- get_min_expr(ligand, counts_matrix)    # ligand

    # bivariate local moran
    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    res_bv[isolate_idx, "Pr(folded) Sim"] <- 1
    res_bv[isolate_idx, "Ibvi"] <- 0

    hs <- spdep::hotspot(
      res_bv,
      Prname = "Pr(folded) Sim",
      cutoff = p_cutoff,
      quadrant.type = "pysal",
      p.adjust = "none"
    )

    idx <- (hs == "High-High")
    idx[is.na(idx)] <- FALSE

    HH_idx  <- which(idx)
    HH_pval <- res_bv[idx, "Pr(folded) Sim"]

    # write back into LR_out
    LR_out$sig_numbers[i] <- length(HH_idx)
    LR_out$sig_index[[i]] <- HH_idx
    LR_out$sig_pval[[i]]  <- HH_pval
    LR_out$all_pval[[i]]  <- res_bv[, "Pr(folded) Sim"]
    LR_out$all_lisa[[i]]  <- res_bv[, "Ibvi"]
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]

  ## ---------------------------
  ## Return both results + binned data
  ## ---------------------------
  front_cols <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")

  LR_out <- LR_out[, c(front_cols, setdiff(colnames(LR_out), front_cols))]

  list(
    LR_out = LR_out,
    bin_sf = bin_sf
  )
}
