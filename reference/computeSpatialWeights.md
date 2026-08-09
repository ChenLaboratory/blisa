# Compute Spatial Weights for BLISA

Builds queen (nearby) and distance-decay (diffuse) spatial weight
matrices from a bin-level `sf` object, excluding bins that have no
neighbours.

## Usage

``` r
computeSpatialWeights(bins, bin_size = 50, dmax = 250)
```

## Arguments

- bins:

  An `sf` object of spatial bins.

- bin_size:

  Numeric. Bin spacing used to define queen adjacency (`1.2 * bin_size`
  radius).

- dmax:

  Numeric. Maximum distance for diffuse-mode neighbours.

## Value

A list with:

- queen_wt:

  Spatial weights list for nearby (queen) mode.

- dist_wt:

  Spatial weights list for diffuse (distance-decay) mode.

- keep_idx_queen:

  Integer indices of bins used in queen-mode Moran.

- keep_idx_dist:

  Integer indices of bins used in diffuse-mode Moran.

- isolate_idx_queen:

  Integer indices of original queen-mode isolates.

- isolate_idx_dist:

  Integer indices of original diffuse-mode isolates.

- queen_nb_full:

  Full (unsubset) neighbour list for nearby mode, indexed over all bins.

- dist_nb_full:

  Full (unsubset) neighbour list for diffuse mode, indexed over all
  bins.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(42)
pts  <- sf::st_as_sf(
  data.frame(x = runif(300, 0, 1000), y = runif(300, 0, 1000)),
  coords = c("x", "y"), crs = NA
)
bins <- sf::st_sf(
  geometry = sf::st_make_grid(pts, cellsize = 100,
                              what = "polygons", square = FALSE)
)
sw <- computeSpatialWeights(bins, bin_size = 100, dmax = 300)
names(sw)
} # }
```
