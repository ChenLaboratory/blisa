hex_binning_cells <- function(
    coords_df = NULL,
    counts_matrix = NULL,
    spe = NULL,
    hex_size = 50) {

  ## ---------------------------
  ## 0. Input handling
  ## ---------------------------
  if (!is.null(spe)) {
    coords_df <- as.data.frame(spe@int_colData$spatialCoords)
    counts_matrix <- spe@assays@data$counts
  } else {
    stopifnot(!is.null(coords_df), !is.null(counts_matrix))
    stopifnot(inherits(counts_matrix, c("matrix", "dgCMatrix")))
  }

  coords_df$cell_id <- rownames(coords_df)

  required_cols <- c("cell_id", "x_centroid", "y_centroid")
  stopifnot(all(required_cols %in% colnames(coords_df)))

  stopifnot(all(colnames(counts_matrix) %in% coords_df$cell_id))

  ## ---------------------------
  ## 1. sf cells
  ## ---------------------------
  cell_sf <- sf::st_as_sf(
    coords_df,
    coords = c("x_centroid", "y_centroid"),
    crs = NA
  )

  ## ---------------------------
  ## 2. Hex binning
  ## ---------------------------
  # create bins
  hex_geom <- sf::st_make_grid( # only polygon info
    cell_sf,
    cellsize = hex_size,
    what = "polygons",
    square = FALSE
  )

  hex_sf <- sf::st_sf( # sf object
    hex_id = seq_along(hex_geom),
    geometry = hex_geom
  )
  dim(hex_sf)

  # map cells to bins
  cell_hex_sf <- sf::st_join(cell_sf, hex_sf, join = sf::st_intersects)
  cell_hex_df <- sf::st_drop_geometry(cell_hex_sf)
  dim(cell_hex_df)
  length(unique(cell_hex_df$hex_id))

  ## ---------------------------
  ## 3. Aggregate counts
  ## ---------------------------
  cell_to_hex <- cell_hex_df$hex_id
  names(cell_to_hex) <- cell_hex_df$cell_id
  cell_to_hex <- cell_to_hex[colnames(counts_matrix)]

  n_hex <- nrow(hex_sf) # hex_sf contains ALL hexes including empty bins
  hex_factor <- factor(cell_to_hex, levels = seq_len(n_hex))

  H <- Matrix::sparse.model.matrix(~ hex_factor - 1)
  hex_gene_counts <- counts_matrix %*% H
  dim(hex_gene_counts)

  hex_sf$n_cells <- as.numeric(Matrix::colSums(H))
  dim(hex_sf)

  list(
    hex_sf = hex_sf,
    hex_gene_counts = hex_gene_counts
  )

}


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
    col, default_mode,
    diffuse_category
  )

  return(res)
}
