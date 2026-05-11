#' Bin cells into hexagonal spatial bins
#'
#' Aggregates single-cell spatial data into hexagonal bins and returns a
#' bin-level count matrix together with a matching \code{sf} polygon object,
#' ready to pass directly to \code{\link{runBLISA.default}}.
#'
#' @param coords_df Data frame or matrix with columns \code{x_centroid} and
#'   \code{y_centroid} (e.g. the output of \code{SpatialExperiment::spatialCoords()}).
#'   Row names must be cell IDs matching the column names of \code{counts_matrix}.
#' @param counts_matrix Gene-by-cell count matrix (dense or sparse). Row names
#'   must be gene symbols; column names must be cell IDs present in
#'   \code{coords_df}.
#' @param bin_size Numeric. Approximate width of each hexagonal bin in
#'   coordinate units (e.g. microns). Analogous to \code{grid.length.x} in
#'   \code{sciderHex::gridDensity}. Default \code{50}.
#' @param min_cells Integer. Bins containing fewer than \code{min_cells} cells
#'   are dropped from the output. Default \code{1}.
#' @param group Factor or character vector of length \code{ncol(counts_matrix)}
#'   giving the cell-type label of each cell. When supplied, a named list of
#'   per-cell-type gene-by-bin matrices is included in the output as
#'   \code{counts_by_group}. Default \code{NULL} (not computed).
#'
#' @return A list with:
#' \describe{
#'   \item{counts_matrix}{Gene-by-bin sparse count matrix (all cells combined).
#'     Column \emph{i} corresponds to row \emph{i} of \code{bin_sf}.}
#'   \item{bin_sf}{An \code{sf} object of hexagonal bin polygons with an
#'     \code{n_cells} column recording how many cells fall in each bin. Row
#'     order matches the columns of \code{counts_matrix}.}
#'   \item{counts_by_group}{(Only when \code{group} is supplied.) A named list
#'     of gene-by-bin sparse matrices, one per cell-type level, with the same
#'     bin order as \code{counts_matrix}.}
#' }
#' @export
hexBinCells <- function(coords_df, counts_matrix, bin_size = 50, min_cells = 1,
                        group = NULL) {

  if (!all(c("x_centroid", "y_centroid") %in% colnames(coords_df)))
    stop("coords_df must contain columns 'x_centroid' and 'y_centroid'.")
  if (!all(colnames(counts_matrix) %in% rownames(coords_df)))
    stop("All column names of counts_matrix must appear as row names in coords_df.")
  if (!is.null(group)) {
    if (length(group) != ncol(counts_matrix))
      stop("'group' must have the same length as ncol(counts_matrix).")
    group <- as.factor(group)
  }

  # Build sf point layer in the same order as counts_matrix columns
  coords_df <- as.data.frame(coords_df)
  cell_sf <- sf::st_as_sf(
    coords_df[colnames(counts_matrix), , drop = FALSE],
    coords = c("x_centroid", "y_centroid"),
    crs    = NA
  )

  # Hexagonal grid covering the full extent of the cells
  hex_geom <- sf::st_make_grid(
    cell_sf,
    cellsize = bin_size,
    what     = "polygons",
    square   = FALSE
  )
  bin_sf <- sf::st_sf(bin_id = seq_along(hex_geom), geometry = hex_geom)

  # Point-in-polygon: assign each cell to its bin (same row order as cell_sf)
  cell_bin_sf <- sf::st_join(cell_sf, bin_sf, join = sf::st_intersects)
  cell_to_bin <- sf::st_drop_geometry(cell_bin_sf)$bin_id

  # Aggregate counts: genes x bins  (H is cells x bins)
  n_bins     <- nrow(bin_sf)
  bin_factor <- factor(cell_to_bin, levels = seq_len(n_bins))
  H          <- Matrix::sparse.model.matrix(~ bin_factor - 1)
  bin_counts <- counts_matrix %*% H

  bin_sf$n_cells <- as.numeric(Matrix::colSums(H))

  # Drop bins below min_cells
  keep       <- bin_sf$n_cells >= min_cells
  bin_sf     <- bin_sf[keep, , drop = FALSE]
  bin_counts <- bin_counts[, keep, drop = FALSE]

  result <- list(counts_matrix = bin_counts, bin_sf = bin_sf)

  # Per-cell-type bin matrices (reuse H, subset rows by cell type)
  if (!is.null(group)) {
    result$counts_by_group <- lapply(levels(group), function(ct) {
      idx <- which(group == ct)
      counts_matrix[, idx, drop = FALSE] %*% H[idx, keep, drop = FALSE]
    })
    names(result$counts_by_group) <- levels(group)
  }

  result
}
