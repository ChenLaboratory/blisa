# Build spot-level hexagonal bins from a Visium SpatialExperiment

Turns a 10x Visium `SpatialExperiment` (as produced by
`scider::readVisium()`) into the pre-binned inputs expected by
[`blisa.default`](https://chenlaboratory.github.io/blisa/reference/blisa.md),
treating each Visium spot as one bin. The spot coordinates are rebuilt
from the integer `array_col`/`array_row` indices onto an exact,
de-tilted hexagonal lattice with 100 um (`spot_pitch`) spacing – the
standalone equivalent of `scider::realignVisium()`.

## Usage

``` r
visiumSpotBins(
  spe,
  spot_pitch = NULL,
  in_tissue_only = TRUE,
  hexagons = TRUE,
  min_total_counts = 10,
  verbose = FALSE
)
```

## Arguments

- spe:

  A Visium `SpatialExperiment` with a `"counts"` assay whose row names
  are gene symbols. `array_col`/`array_row` in `colData` are used to
  rebuild the exact lattice when present; otherwise `spatialCoords(spe)`
  is used as-is.

- spot_pitch:

  Numeric or `NULL`. Spot center-to-center spacing used as `bin_size`.
  `NULL` (default) means: use 100 um for the array-index lattice, or
  measure it from the coordinates in the fallback. Supply a number to
  force it (interpreted in coordinate units).

- in_tissue_only:

  Logical. Keep only spots with `in_tissue == 1`. Default `TRUE`.

- hexagons:

  Logical. If `TRUE` (default), build hexagonal bin polygons; if
  `FALSE`, use point geometry (faster, but `plotHotspots` then draws
  points).

- min_total_counts:

  Numeric. Drop spots whose total counts fall below this threshold.
  Default `10`. Set to `0` to keep all spots.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

## Value

A list with:

- counts_matrix:

  Gene-by-spot count matrix (columns match `bins`).

- bins:

  An `sf` object of per-spot bins with `bin_id`, `n_cells` (always 1),
  `total_counts`, and `img_x`/ `img_y` (the raw image-registered spot
  coordinates, for overlaying plots on the H&E image) columns. Row order
  matches the columns of `counts_matrix`.

- pitch:

  The spot spacing used (= the recommended `bin_size`).

- coord_unit:

  `"micron"`, `"pixel/other"`, or `NA` (when `spot_pitch` was supplied)
  – the inferred coordinate unit.

## Details

On Visium, each spot is already a spatial bin, so no cell-to-bin
aggregation is needed (unlike
[`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)).
Because `array_col` increments by 2 within a row and alternate rows
shift `array_col` by 1, the half-row stagger is already encoded in the
indices; the lattice therefore needs no odd/even phase correction. Only
the bin centroids are used by
[`computeSpatialWeights`](https://chenlaboratory.github.io/blisa/reference/computeSpatialWeights.md);
the hexagon polygons are for plotting (they do tessellate the lattice
exactly).

When `array_col`/`array_row` are absent (e.g. a `SpatialExperiment` not
built by `scider::readVisium()`), the function falls back to using
`spatialCoords(spe)` directly: spots are still one bin each, but the
de-tilt/exact-lattice step is skipped. In that case, if `spot_pitch` is
`NULL`, the pitch is measured from the coordinates as the median
nearest-neighbour distance. Since the physical Visium pitch is 100 um,
`100 / pitch` is microns-per-unit; the coordinates are *rescaled to
microns* by that factor (near-identity when they are already microns),
so the returned `pitch` is 100 and downstream `bin_size`/`dmax` are in
microns. If you instead supply `spot_pitch`, no rescaling is done and
everything stays in the coordinates' native units.

## See also

[`blisa`](https://chenlaboratory.github.io/blisa/reference/blisa.md),
[`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
