# Create a blisa object

Constructor for the `blisa` S3 class, which stores the full output of
[`blisa`](https://chenlaboratory.github.io/blisa/reference/blisa.md).

## Usage

``` r
new_blisa(
  LR_results,
  bins,
  spatial_weights,
  CCI_scores = NULL,
  level = c("lr", "pathway")
)
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

- level:

  Character. The unit each row of `LR_results` represents: `"lr"`
  (ligand-receptor pairs, the default) or `"pathway"` (from
  [`blisaPathway`](https://chenlaboratory.github.io/blisa/reference/blisaPathway.md)).
  Controls labelling in `print` and plots.

## Value

An object of class `blisa`.
