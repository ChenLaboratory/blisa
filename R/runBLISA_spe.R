
#' Run BLISA on a cell-level SpatialExperiment
#'
#' Bins cells into hexagonal tiles via \code{hex_binning_cells}, removes empty
#' bins, then runs \code{runBLISA.default}.
#'
#' @param spe A cell-level \code{SpatialExperiment} object with spatial
#'   coordinates in \code{spatialCoords} and counts in \code{assay(spe, "counts")}.
#' @param sf Unused; reserved for future gridding overrides.
#' @param LR_df Data frame of ligand-receptor pairs (see \code{filterLRpairs}).
#' @param hex_size Numeric. Hexagonal bin side length in coordinate units.
#'   Default 50.
#' @param dmax Numeric. Maximum distance for diffuse-mode neighbours. Default 250.
#' @param nsim Integer. Number of permutations for Moran's I significance.
#'   Default 999.
#' @param p_cutoff Numeric. P-value threshold for High-High hotspots. Default 0.05.
#' @param min_ligand Numeric. Minimum ligand count threshold. Default 10.
#' @param min_receptor Numeric. Minimum receptor count threshold. Default 10.
#' @param col Character. Column in \code{LR_df} for interaction category.
#'   Default \code{"annotation"}.
#' @param default_mode Character. Default CCC mode when \code{col} is absent.
#'   Default \code{"diffuse"}.
#' @param diffuse_category Character vector of annotation categories treated as
#'   diffuse signaling.
#'
#' @return A list with \code{LR_out} and \code{bin_sf}; see
#'   \code{\link{runBLISA.default}} for details.
#' @export
runBLISA.spe <- function(
    spe, # default cell -> bining inside
    sf = NULL, # gridding info
    LR_df, # default cellchatDB (download to local/ git online read)
    hex_size = 50,
    dmax = 250,
    nsim = 999,
    p_cutoff = 0.05,
    min_ligand = 10, min_receptor = 10,
    col = "annotation", default_mode = "diffuse",
    diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
) {
  coords_df <- as.data.frame(spe@int_colData$spatialCoords)
  counts_matrix <- spe@assays@data$counts

  ## ---------------------------
  ## Hex binning
  ## ---------------------------
  binning_res <- hex_binning_cells(coords_df, counts_matrix, spe, hex_size)
  hex_gene_counts <- binning_res$hex_gene_counts
  hex_sf <- binning_res$hex_sf

  # Remove empty bins
  keep_idx <- which(hex_sf$n_cells > 0)

  hex_sf <- hex_sf[keep_idx, ]
  hex_gene_counts <- hex_gene_counts[, keep_idx, drop = FALSE]

  ## ---------------------------
  ## run BLISA
  ## ---------------------------

  res <- runBLISA.default( # bin-level
    counts_matrix = hex_gene_counts,
    bin_sf = hex_sf,
    LR_df, # use index, remove other info
    hex_size,
    dmax,
    nsim,
    p_cutoff,
    min_ligand, min_receptor,
    col, default_mode,
    diffuse_category
  )

  return(res)
}

#' Run BLISA on a cell-level SpatialExperiment with isolate removal
#'
#' Bins cells into hexagonal tiles via \code{hex_binning_cells}, removes empty
#' bins, then runs \code{runBLISA.default.isolates.removed}. Isolated bins and
#' bins below \code{min_cells_per_bin} are excluded from Moran's I and assigned
#' neutral statistics (p = 1, LISA = 0).
#'
#' @param spe A cell-level \code{SpatialExperiment} object.
#' @param sf Unused; reserved for future gridding overrides.
#' @param LR_df Data frame of ligand-receptor pairs (see \code{filterLRpairs}).
#' @param hex_size Numeric. Hexagonal bin side length. Default 50.
#' @param dmax Numeric. Maximum distance for diffuse-mode neighbours. Default 250.
#' @param nsim Integer. Number of permutations for Moran's I. Default 999.
#' @param p_cutoff Numeric. P-value threshold for hotspots. Default 0.05.
#' @param min_ligand Numeric. Minimum ligand count threshold. Default 10.
#' @param min_receptor Numeric. Minimum receptor count threshold. Default 10.
#' @param min_cells_per_bin Integer. Minimum number of cells a hex bin must
#'   contain to be included in spatial weight computation and Moran's I
#'   calculation. Bins below this threshold are excluded and assigned neutral
#'   statistics (p = 1, LISA = 0). Default 1.
#' @param col Character. Column in \code{LR_df} for interaction category.
#'   Default \code{"annotation"}.
#' @param default_mode Character. Default CCC mode when \code{col} is absent.
#'   Default \code{"diffuse"}.
#' @param diffuse_category Character vector of annotation categories treated as
#'   diffuse signaling.
#'
#' @return A list with \code{LR_out} and \code{bin_sf}; see
#'   \code{\link{runBLISA.default.isolates.removed}} for details.
#' @export
runBLISA.spe.isolates.removed <- function(
    spe, # default cell -> bining inside
    sf = NULL, # gridding info
    LR_df, # default cellchatDB (download to local/ git online read)
    hex_size = 50,
    dmax = 250,
    nsim = 999,
    p_cutoff = 0.05,
    min_ligand = 10, min_receptor = 10,
    min_cells_per_bin = 1,
    col = "annotation", default_mode = "diffuse",
    diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
) {
  coords_df <- as.data.frame(spe@int_colData$spatialCoords)
  counts_matrix <- spe@assays@data$counts

  ## ---------------------------
  ## Hex binning
  ## ---------------------------
  binning_res <- hex_binning_cells(coords_df, counts_matrix, spe, hex_size)
  hex_gene_counts <- binning_res$hex_gene_counts
  hex_sf <- binning_res$hex_sf

  # Remove empty bins
  keep_idx <- which(hex_sf$n_cells > 0)

  hex_sf <- hex_sf[keep_idx, ]
  hex_gene_counts <- hex_gene_counts[, keep_idx, drop = FALSE]

  ## ---------------------------
  ## run BLISA
  ## ---------------------------

  res <- runBLISA.default.isolates.removed( # bin-level
    counts_matrix = hex_gene_counts,
    bin_sf = hex_sf,
    LR_df, # use index, remove other info
    hex_size,
    dmax,
    nsim,
    p_cutoff,
    min_ligand, min_receptor,
    min_cells_per_bin,
    n_cells_col = "n_cells",
    col, default_mode,
    diffuse_category
  )

  return(res)
}
