
#' Run BLISA analysis on cell-level spe
#'
#' @param spe cell-level spe object
#' @param bin_sf Bin-level sf object
#' @param LR_df Ligand-receptor dataframe
#' @export
#'
#'
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

#' Run BLISA on Bin-Level Spatial Data with Isolate Removal
#'
#' @param min_cells_per_bin Integer. Minimum number of cells a hex bin must
#'   contain to be included in spatial weight computation and Moran's I
#'   calculation. Bins below this threshold are excluded and assigned neutral
#'   statistics (p = 1, LISA = 0).
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
