# Sender-by-receiver heatmap of CCI scores for one ligand-receptor pair

Generic function. Reshapes the CCI data frame into a receiver-by-sender
cell-type matrix for one selected LR pair and draws a clustered heatmap.

## Usage

``` r
plotCCILR(x, ...)

# S3 method for class 'blisa'
plotCCILR(x, index = 1, ligand = NULL, receptor = NULL, main = NULL, ...)

# Default S3 method
plotCCILR(x, lr_pair, main = NULL, ...)
```

## Arguments

- x:

  A `blisa` object or a CCI scores data frame (the `CCI_scores` slot of
  a `blisa` object).

- ...:

  Additional arguments passed to the relevant method.

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

- main:

  Character or `NULL`. Title drawn above the heatmap. When supplied, the
  heatmap is drawn with this overall title (via
  [`ComplexHeatmap::draw`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html));
  the `Heatmap` object is returned invisibly. Default `NULL` (no title;
  object returned for the caller to print).

- lr_pair:

  Character. Column name in the CCI scores data frame corresponding to
  the ligand-receptor pair to visualise (e.g. `"CXCL12_CXCR4"`).

## Value

A `Heatmap` object.

## Methods (by class)

- `plotCCILR(blisa)`: Method for a `blisa` object. The LR pair is
  selected by `index` (default 1, the top-ranked pair) unless both
  `ligand` and `receptor` are supplied, in which case the matching row
  is located automatically and `index` is ignored. Stops with an
  informative error if `CCI_scores` is `NULL` or the selected LR pair
  has no significant hotspots.

- `plotCCILR(default)`: Method for a CCI scores data frame (e.g. the
  `CCI_scores` slot of a `blisa` object). The LR pair is selected by
  column name via `lr_pair`.

## See also

[`plotCCI`](https://chenlaboratory.github.io/blisa/reference/plotCCI.md)
for an overview heatmap across all LR pairs;
[`plotCCIsummary`](https://chenlaboratory.github.io/blisa/reference/plotCCIsummary.md)
for an aggregated sender-by-receiver heatmap.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")
plotCCILR(result, index = 1)
} # }
```
