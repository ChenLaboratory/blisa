# Bin cells into hexagonal spatial bins

Aggregates single-cell spatial data into hexagonal bins and returns a
bin-level count matrix together with a matching `sf` polygon object,
ready to pass directly to
[`blisa.default`](https://chenlaboratory.github.io/blisa/reference/blisa.md).

## Usage

``` r
hexBinCells(
  coords_df,
  counts_matrix,
  bin_size = 50,
  min_cells = 1,
  min_total_counts = 10,
  group = NULL,
  verbose = FALSE
)
```

## Arguments

- coords_df:

  Data frame or matrix with columns `x_centroid` and `y_centroid` (e.g.
  the output of
  [`SpatialExperiment::spatialCoords()`](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment-methods.html)).
  Row names must be cell IDs matching the column names of
  `counts_matrix`.

- counts_matrix:

  Gene-by-cell count matrix (dense or sparse). Row names must be gene
  symbols; column names must be cell IDs present in `coords_df`.

- bin_size:

  Numeric. Approximate width of each hexagonal bin in coordinate units
  (e.g. microns). Analogous to `grid.length.x` in
  `sciderHex::gridDensity`. Default `50`.

- min_cells:

  Integer. Bins containing fewer than `min_cells` cells are dropped from
  the output. Default `1`.

- min_total_counts:

  Numeric. Bins whose total counts (summed over all genes) fall below
  this threshold are dropped from the output, alongside the `min_cells`
  filter. Set to `0` to disable. Default `10`.

- group:

  Factor or character vector of length `ncol(counts_matrix)` giving the
  cell-type label of each cell. When supplied, a named list of
  per-cell-type gene-by-bin matrices is included in the output as
  `counts_by_group`. Default `NULL` (not computed).

- verbose:

  Logical. If `TRUE`, print progress messages. Default `FALSE` (silent).

## Value

A list with:

- counts_matrix:

  Gene-by-bin sparse count matrix (all cells combined). Column *i*
  corresponds to row *i* of `bins`.

- bins:

  An `sf` object of hexagonal bin polygons with an `n_cells` column
  recording how many cells fall in each bin and a `total_counts` column
  recording the summed counts per bin. Row order matches the columns of
  `counts_matrix`.

- counts_by_group:

  (Only when `group` is supplied.) A named list of gene-by-bin sparse
  matrices, one per cell-type level, with the same bin order as
  `counts_matrix`.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
n <- 500
coords <- data.frame(
  x_centroid = runif(n, 0, 1000),
  y_centroid = runif(n, 0, 1000),
  row.names  = paste0("cell_", seq_len(n))
)
counts <- Matrix::Matrix(
  matrix(rpois(20L * n, lambda = 5), nrow = 20L, ncol = n,
         dimnames = list(paste0("gene_", 1:20), paste0("cell_", seq_len(n)))),
  sparse = TRUE
)
group  <- sample(c("TypeA", "TypeB"), n, replace = TRUE)
binned <- hexBinCells(coords, counts, bin_size = 100, group = group)
str(binned, max.level = 1)
} # }
```
