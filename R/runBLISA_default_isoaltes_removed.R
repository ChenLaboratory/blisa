#' Run BLISA on Bin-Level Spatial Data with Isolate Removal
#'
#' Performs BLISA (Bivariate Local Indicator of Spatial Association) analysis
#' on bin-level spatial transcriptomics data to identify spatially enriched
#' ligand-receptor interactions. Spatially isolated bins and bins failing
#' minimum cell count filtering are excluded from Moran's I calculation and
#' assigned neutral statistics in the final output (p = 1, I = 0).
#'
#' @param counts_matrix Numeric gene-by-bin expression matrix. Rows correspond
#'   to genes and columns correspond to spatial bins.
#' @param bin_sf An \code{sf} object containing bin-level spatial geometries
#'   and metadata. Must contain one row per bin.
#' @param LR_df Data frame of ligand-receptor pairs. Must contain ligand and
#'   receptor gene symbol columns compatible with \code{filterLRpairs()}.
#' @param hex_size Numeric. Approximate bin/hexagon spacing used to define
#'   queen adjacency for nearby signaling mode.
#' @param dmax Numeric. Maximum diffusion distance (in spatial coordinate units)
#'   used for diffuse signaling mode.
#' @param nsim Integer. Number of permutations for local bivariate Moran's I
#'   significance estimation.
#' @param p_cutoff Numeric. P-value threshold for calling significant
#'   High-High ligand-receptor hotspots.
#' @param min_ligand Numeric. Minimum ligand expression threshold required for
#'   a ligand to be retained in ligand-receptor filtering.
#' @param min_receptor Numeric. Minimum receptor expression threshold required
#'   for a receptor to be retained in ligand-receptor filtering.
#' @param min_cells_per_bin Integer. Minimum number of cells required for a bin
#'   to be included in Moran's I calculation. Bins below this threshold are
#'   excluded and assigned neutral statistics.
#' @param col Character. Column name in \code{LR_df} specifying ligand-receptor
#'   annotation/category used for communication mode assignment.
#' @param default_mode Character. Default communication mode assigned to ligand-
#'   receptor pairs not matching \code{diffuse_category}. Typically one of
#'   \code{"nearby"} or \code{"diffuse"}.
#' @param diffuse_category Character vector of ligand-receptor annotation
#'   categories to be treated as diffuse signaling.
#'
#' @returns A list containing:
#' \describe{
#'   \item{LR_out}{Data frame summarizing BLISA results for each ligand-receptor
#'   pair, including significant hotspot counts, hotspot indices, p-values,
#'   and local Moran's I statistics.}
#'   \item{bin_sf}{Input bin-level \code{sf} object.}
#'   \item{isolate_idx_queen}{Integer vector of bins with no queen-adjacent
#'   neighbors.}
#'   \item{isolate_idx_dist}{Integer vector of bins with no neighbors within
#'   \code{dmax}.}
#'   \item{low_cell_idx}{Integer vector of bins excluded for having fewer than
#'   \code{min_cells_per_bin} cells.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- runBLISA.default.isoaltes.removed(
#'   counts_matrix = counts_mat,
#'   bin_sf = bin_sf,
#'   LR_df = LR_pairs
#' )
#' }
runBLISA.default.isolates.removed <- function(
    counts_matrix,
    bin_sf,
    LR_df,
    hex_size = 50,
    dmax = 250,
    nsim = 999,
    p_cutoff = 0.05,
    min_ligand = 10, min_receptor = 10,
    min_cells_per_bin = 1,
    col = "annotation",
    default_mode = "diffuse",
    diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
) {

  centroids <- sf::st_centroid(bin_sf)
  coords <- sf::st_coordinates(centroids)

  ## ---------------------------
  ## Filter low-cell bins
  ## ---------------------------

  if (!"n_cells" %in% colnames(bin_sf)) {
    stop("bin_sf must contain a 'n_cells' column for min_cells_per_bin filtering.")
  }

  low_cell_idx <- which(bin_sf$n_cells < min_cells_per_bin)

  message(length(low_cell_idx),
          " bins removed for having < ",
          min_cells_per_bin,
          " cells.")

  ## ---------------------------
  ## Queen spatial weights
  ## ---------------------------

  # queen neighboours of all bins
  queen_nb_full <- spdep::dnearneigh(coords, 0, 1.2 * hex_size)

  # isolates  bins with no queen neighboours
  isolate_idx_queen <- which(spdep::card(queen_nb_full) == 0)

  message(length(isolate_idx_queen),
          " isolated bins with no nearby neighbors: ",
          paste(isolate_idx_queen, collapse = ","))

  # idx to use for queen weights: exclude isolates and low-cell bins
  exclude_idx_queen <- union(isolate_idx_queen, low_cell_idx)
  keep_idx_queen <- setdiff(seq_len(nrow(bin_sf)), exclude_idx_queen)
  keep_logical_queen <- seq_len(nrow(bin_sf)) %in% keep_idx_queen

  # queen neighboours of kept  bins
  queen_nb <- spdep::subset.nb(
    queen_nb_full,
    subset = keep_logical_queen
  )

  # queen weights of kept bins
  queen_wt <- spdep::nb2listwdist(
    queen_nb,
    centroids[keep_idx_queen, ],
    type = "idw",
    style = "W",
    zero.policy = TRUE
  )

  ## ---------------------------
  ## Distance spatial weights
  ## ---------------------------

  # distance neighboours of all bins
  dist_nb_full <- spdep::dnearneigh(coords, 0, dmax)

  # isolates  bins with no distance neighboours
  isolate_idx_dist <- which(spdep::card(dist_nb_full) == 0)
  message(length(isolate_idx_dist),
          " isolated bins with no neighbors within ", dmax, " um: ",
          paste(isolate_idx_dist, collapse = ","))

  # idx to use for distance weights: exclude isolates and low-cell bins
  exclude_idx_dist  <- union(isolate_idx_dist, low_cell_idx)

  keep_idx_dist  <- setdiff(seq_len(nrow(bin_sf)), exclude_idx_dist)
  keep_logical_dist  <- seq_len(nrow(bin_sf)) %in% keep_idx_dist

  # distance neighboours of non-isolates  bins
  dist_nb <- spdep::subset.nb(dist_nb_full, subset = keep_logical_dist)

  # distance weights of non-isolates  bins
  weight_at_dmax <- 0.01

  dist_wt <- spdep::nb2listwdist(
    dist_nb,
    bin_sf[keep_idx_dist, ],
    type = "exp",
    style = "W",
    alpha = -log(weight_at_dmax) / dmax
  )

  ## ---------------------------
  ## Filter LR pairs
  ## ---------------------------
  LR_df_filtered <- filterLRpairs(
    counts = counts_matrix,
    min_ligand,
    min_receptor,
    LR_df
  )

  ## ---------------------------
  ## Local bivariate Moran
  ## ---------------------------
  LR_out <- LR_df_add_mode(
    LR_df_filtered,
    col,
    default_mode,
    diffuse_category
  )

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

    if (mode == "nearby") {
      wt <- queen_wt
      keep_idx <- keep_idx_queen
      isolate_idx <- isolate_idx_queen
    } else {
      wt <- dist_wt
      keep_idx <- keep_idx_dist
      isolate_idx <- isolate_idx_dist
    }

    message("ccc mode is ", mode)

    x_full <- get_min_expr(receptor, counts_matrix)
    y_full <- get_min_expr(ligand, counts_matrix)

    x <- x_full[keep_idx]
    y <- y_full[keep_idx]

    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    ## Restore full-length vectors
    full_pval <- rep(1, length(x_full))
    full_lisa <- rep(0, length(x_full))

    full_pval[keep_idx] <- res_bv[, "Pr(folded) Sim"]
    full_lisa[keep_idx] <- res_bv[, "Ibvi"]

    hs <- spdep::hotspot(
      res_bv,
      Prname = "Pr(folded) Sim",
      cutoff = p_cutoff,
      quadrant.type = "pysal",
      p.adjust = "none"
    )

    idx_keep <- (hs == "High-High")
    idx_keep[is.na(idx_keep)] <- FALSE

    HH_idx <- keep_idx[which(idx_keep)]
    HH_pval <- full_pval[HH_idx]

    LR_out$sig_numbers[i] <- length(HH_idx)
    LR_out$sig_index[[i]] <- HH_idx
    LR_out$sig_pval[[i]]  <- HH_pval
    LR_out$all_pval[[i]]  <- full_pval
    LR_out$all_lisa[[i]]  <- full_lisa
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]

  front_cols <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")

  LR_out <- LR_out[, c(front_cols, setdiff(colnames(LR_out), front_cols))]

  list(
    LR_out = LR_out,
    bin_sf = bin_sf,
    isolate_idx_queen = isolate_idx_queen,
    isolate_idx_dist = isolate_idx_dist
  )
}
