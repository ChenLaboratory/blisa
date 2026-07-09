#' Bin cells into hexagonal spatial bins
#'
#' Aggregates single-cell spatial data into hexagonal bins and returns a
#' bin-level count matrix together with a matching \code{sf} polygon object,
#' ready to pass directly to \code{\link{blisa.default}}.
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
#' @param min_total_counts Numeric. Bins whose total counts (summed over all
#'   genes) fall below this threshold are dropped from the output, alongside
#'   the \code{min_cells} filter. Set to \code{0} to disable. Default \code{10}.
#' @param group Factor or character vector of length \code{ncol(counts_matrix)}
#'   giving the cell-type label of each cell. When supplied, a named list of
#'   per-cell-type gene-by-bin matrices is included in the output as
#'   \code{counts_by_group}. Default \code{NULL} (not computed).
#'
#' @return A list with:
#' \describe{
#'   \item{counts_matrix}{Gene-by-bin sparse count matrix (all cells combined).
#'     Column \emph{i} corresponds to row \emph{i} of \code{bins}.}
#'   \item{bins}{An \code{sf} object of hexagonal bin polygons with an
#'     \code{n_cells} column recording how many cells fall in each bin and a
#'     \code{total_counts} column recording the summed counts per bin. Row
#'     order matches the columns of \code{counts_matrix}.}
#'   \item{counts_by_group}{(Only when \code{group} is supplied.) A named list
#'     of gene-by-bin sparse matrices, one per cell-type level, with the same
#'     bin order as \code{counts_matrix}.}
#' }
#' @examples
#' \dontrun{
#' set.seed(42)
#' n <- 500
#' coords <- data.frame(
#'   x_centroid = runif(n, 0, 1000),
#'   y_centroid = runif(n, 0, 1000),
#'   row.names  = paste0("cell_", seq_len(n))
#' )
#' counts <- Matrix::Matrix(
#'   matrix(rpois(20L * n, lambda = 5), nrow = 20L, ncol = n,
#'          dimnames = list(paste0("gene_", 1:20), paste0("cell_", seq_len(n)))),
#'   sparse = TRUE
#' )
#' group  <- sample(c("TypeA", "TypeB"), n, replace = TRUE)
#' binned <- hexBinCells(coords, counts, bin_size = 100, group = group)
#' str(binned, max.level = 1)
#' }
#' @export
hexBinCells <- function(coords_df, counts_matrix, bin_size = 50, min_cells = 1,
                        min_total_counts = 10, group = NULL) {

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
  bins <- sf::st_sf(bin_id = seq_along(hex_geom), geometry = hex_geom)

  # Point-in-polygon: assign each cell to its bin (same row order as cell_sf).
  # A cell whose coordinates fall exactly on a shared hexagon boundary can
  # intersect two bins; st_join then emits two rows for that cell, making H
  # non-conformable with counts_matrix. Stamp each cell with a stable index
  # before joining and deduplicate afterward, keeping the first bin match.
  cell_sf$.cell_idx <- seq_len(nrow(cell_sf))
  cell_bins    <- sf::st_join(cell_sf, bins, join = sf::st_intersects)
  cell_bins_df <- sf::st_drop_geometry(cell_bins)
  cell_bins_df <- cell_bins_df[!duplicated(cell_bins_df$.cell_idx), , drop = FALSE]
  cell_to_bin  <- cell_bins_df$bin_id
  cell_sf$.cell_idx <- NULL

  # Aggregate counts: genes x bins  (H is cells x bins)
  n_bins     <- nrow(bins)
  bin_factor <- factor(cell_to_bin, levels = seq_len(n_bins))
  H          <- Matrix::sparse.model.matrix(~ bin_factor - 1)
  bin_counts <- counts_matrix %*% H

  bins$n_cells      <- as.numeric(Matrix::colSums(H))
  bins$total_counts <- as.numeric(Matrix::colSums(bin_counts))

  # Drop bins below min_cells or below min_total_counts
  keep       <- bins$n_cells >= min_cells & bins$total_counts >= min_total_counts
  message(sum(!keep), " bins dropped (< ", min_cells, " cells or < ",
          min_total_counts, " total counts).")
  bins       <- bins[keep, , drop = FALSE]
  bin_counts <- bin_counts[, keep, drop = FALSE]

  result <- list(counts_matrix = bin_counts, bins = bins)

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
