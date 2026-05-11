#' Bin Single Cells into Hexagonal Tiles (Deprecated)
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is superseded by \code{\link{hexBinCells}}, which has a
#' cleaner interface, supports per-group binning via the \code{group} argument,
#' and returns output compatible with \code{\link{runBLISA}}.
#'
#' @param coords_df Data frame with columns \code{x_centroid}, \code{y_centroid},
#'   and row names as cell IDs. Required when \code{spe} is \code{NULL}.
#' @param counts_matrix Gene-by-cell count matrix (dense or sparse). Required
#'   when \code{spe} is \code{NULL}.
#' @param spe A \code{SpatialExperiment} object. When provided, \code{coords_df}
#'   and \code{counts_matrix} are extracted automatically.
#' @param hex_size Numeric. Side length of each hexagonal bin in coordinate units.
#'
#' @return A list with \code{hex_sf}, \code{hex_gene_counts}, and
#'   \code{cell_to_hex}. Use \code{\link{hexBinCells}} for new code.
#' @export
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

  ## ---------------------------
  ## 3. Aggregate counts
  ## ---------------------------
  cell_to_hex <- setNames(cell_hex_df$hex_id, cell_hex_df$cell_id)
  cell_to_hex <- cell_to_hex[colnames(counts_matrix)]

  n_hex <- nrow(hex_sf) # hex_sf contains ALL hexes including empty bins
  hex_factor <- factor(cell_to_hex, levels = seq_len(n_hex))
  H <- Matrix::sparse.model.matrix(~ hex_factor - 1)

  hex_gene_counts <- counts_matrix %*% H
  hex_sf$n_cells <- as.numeric(Matrix::colSums(H))

  list(
    hex_sf = hex_sf,
    hex_gene_counts = hex_gene_counts,
    cell_to_hex     = cell_to_hex
  )

}

