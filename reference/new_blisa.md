# Create a blisa object

Constructor for the `blisa` S3 class, which stores the full output of
[`blisa`](https://chenlaboratory.github.io/blisa/reference/blisa.md).

## Usage

``` r
new_blisa(LR_results, bins, spatial_weights, CCI_scores = NULL)
```

## Arguments

- LR_results:

  Data frame of BLISA results, one row per ligand-receptor pair.

- bins:

  An `sf` object of hexagonal bin polygons.

- spatial_weights:

  Spatial weights list from
  [`computeSpatialWeights`](https://chenlaboratory.github.io/blisa/reference/computeSpatialWeights.md).

- CCI_scores:

  Wide data frame of cell-cell interaction scores from
  [`runCCI`](https://chenlaboratory.github.io/blisa/reference/runCCI.md),
  or `NULL` if CCI was not computed.

## Value

An object of class `blisa`.
