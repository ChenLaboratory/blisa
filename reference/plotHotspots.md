# Spatial hotspot map for one ligand-receptor pair

Generic function. Plots each bin coloured by significance status: empty,
non-significant, or significant hotspot (continuous gradient of -log10
p-value or 1 - p-value).

## Usage

``` r
plotHotspots(x, ...)

# S3 method for class 'blisa'
plotHotspots(
  x,
  index = 1,
  ligand = NULL,
  receptor = NULL,
  as_points = FALSE,
  size = 1.5,
  background = NULL,
  log_pval = TRUE,
  p_cutoff = NULL,
  spots = NULL,
  spot_pval = NULL,
  title = NULL,
  ...
)
```

## Arguments

- x:

  A `blisa` object.

- ...:

  Additional arguments passed to the method.

- index:

  Integer. Row index into `LR_results` selecting the ligand-receptor
  pair to visualise. Ignored when both `ligand` and `receptor` are
  supplied. Default `1` (top-ranked pair).

- ligand:

  Character. Ligand gene symbol. When both `ligand` and `receptor` are
  provided the matching LR pair is located automatically and `index` is
  ignored. Must be supplied together with `receptor`.

- receptor:

  Character. Receptor gene symbol. Must be supplied together with
  `ligand`.

- as_points:

  Logical. If `TRUE`, draw each bin as a dot at its centroid instead of
  its polygon. Useful for Visium, where each bin is a single spot.
  Default `FALSE` (draw polygons).

- size:

  Numeric. Point size when `as_points = TRUE` (or when `background` is
  supplied). Default `1.5`.

- background:

  A `ggplot` object to draw the hotspots on top of, e.g.
  `scider::plotImage(spe)` to place the H&E image behind the spots. When
  supplied, the hotspots are always rendered as dots (the raster
  background is incompatible with `geom_sf` coordinates), and the bins
  must share the background plot's coordinate frame (true for objects
  read by `scider::readVisium()`). Default `NULL`.

- log_pval:

  Logical. If `TRUE` (default), colour significant bins by
  -log10(p-value). If `FALSE`, use 1 - p-value.

- p_cutoff:

  Numeric or `NULL`. When `NULL` (default), the pre-computed hotspot
  bins stored in the `blisa` object are used, reflecting the `p_cutoff`
  and High-High quadrant classification applied during
  [`blisa`](https://chenlaboratory.github.io/blisa/reference/blisa.md).
  When a numeric value is supplied, bins are re-defined on the fly as
  those with `all_pval <= p_cutoff` and quadrant label `"High-High"`
  (from the stored `all_quadrant`), giving an exact re-threshold
  consistent with the original classification.

- spots:

  Integer vector or `NULL`. Significant spot indices to plot directly,
  bypassing the LR-pair lookup – e.g. a `sig_spot_index` entry from
  [`blisaPathway`](https://chenlaboratory.github.io/blisa/reference/blisaPathway.md).
  When supplied, `index`/ `ligand`/`receptor`/`p_cutoff` are ignored and
  every bin is drawn as grey background with these spots coloured.
  Default `NULL`.

- spot_pval:

  Numeric vector. P-values for `spots` (same length and order); required
  when `spots` is supplied.

- title:

  Character or `NULL`. Plot title. Defaults to the LR pair id (or
  "hotspots" in `spots` mode). Useful to label a pathway, e.g.
  `title = "MDK pathway"`.

## Value

A `ggplot` object.

## Methods (by class)

- `plotHotspots(blisa)`: Method for a `blisa` object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")
plotHotspots(result, index = 1)
plotHotspots(result, index = 1, log_pval = FALSE)
} # }
```
