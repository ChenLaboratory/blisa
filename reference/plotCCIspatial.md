# Spatial map of dominant sender-receiver cell-type pairs at BLISA hotspots

For a selected ligand-receptor pair, identifies the dominant interacting
cell-type pair at each significant hotspot bin and draws a spatial map
of the tissue coloured by those pairs. Receiver cells are those inside
hotspot bins; sender cells are drawn from the immediate neighbourhood.

## Usage

``` r
plotCCIspatial(
  x,
  counts_by_group,
  index = 1,
  ligand = NULL,
  receptor = NULL,
  top_pairs = 30
)
```

## Arguments

- x:

  A `blisa` object as returned by
  [`blisa`](https://chenlaboratory.github.io/blisa/reference/blisa.md).

- counts_by_group:

  Named list of gene-by-bin count matrices, one per cell type. Typically
  the `counts_by_group` element returned by
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md).
  Names must match the cell-type levels.

- index:

  Integer. Row index into `x$LR_results` selecting the ligand-receptor
  pair to visualise. Ignored when both `ligand` and `receptor` are
  supplied. Default `1` (top-ranked pair).

- ligand:

  Character. Ligand gene symbol. When both `ligand` and `receptor` are
  provided the matching LR pair is located automatically and `index` is
  ignored. Must be supplied together with `receptor`.

- receptor:

  Character. Receptor gene symbol. Must be supplied together with
  `ligand`.

- top_pairs:

  Integer. Maximum number of distinct cell-type pairs to show in the
  legend; remaining pairs are grouped as `"rare pairs"`. Default `30`.

## Value

A `ggplot` object.

## See also

[`plotHotspots`](https://chenlaboratory.github.io/blisa/reference/plotHotspots.md)
for a significance-based spatial map of hotspot bins.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")
binned <- hexBinCells(
  as.data.frame(SpatialExperiment::spatialCoords(spe)),
  SummarizedExperiment::assay(spe, "counts"),
  bin_size = 50, group = spe$cell_type
)
plotCCIspatial(result, binned$counts_by_group, index = 1)
} # }
```
