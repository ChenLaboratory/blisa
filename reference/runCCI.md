# Score cell-cell interactions from BLISA hotspots

Generic function for scoring cell-cell interactions. Dispatches on the
class of `x`:

- `runCCI.blisa` accepts a `blisa` object. If `CCI_scores` are already
  present and `overwrite = FALSE` (the default), the object is returned
  unchanged. Set `overwrite = TRUE` with a `counts_by_group` to
  recompute and replace existing scores. If no scores exist,
  `counts_by_group` must be supplied and scores are computed and
  attached.

- `runCCI.default` performs the raw computation given a `blisa` object
  and a `counts_by_group` list, returning only the scores data frame.
  Used internally by `runCCI.blisa` and
  [`blisa.default`](https://chenlaboratory.github.io/blisa/reference/blisa.md).

## Usage

``` r
runCCI(x, ...)

# S3 method for class 'blisa'
runCCI(x, counts_by_group = NULL, overwrite = FALSE, ...)

# Default S3 method
runCCI(x, counts_by_group, ...)
```

## Arguments

- x:

  A `blisa` object.

- ...:

  Additional arguments passed to the relevant method.

- counts_by_group:

  Named list of gene-by-bin sparse count matrices, one per group level
  (e.g. cell type). Typically the `counts_by_group` element of the list
  returned by
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  when `group` is supplied. Names must match the group levels. Required
  when `x$CCI_scores` is `NULL` or when `overwrite = TRUE`.

- overwrite:

  Logical. If `FALSE` (default) and `x$CCI_scores` is already populated,
  the object is returned unchanged. If `TRUE` and `counts_by_group` is
  supplied, existing scores are recomputed and replaced.

## Value

See individual method documentation.

`runCCI.blisa`: the input `blisa` object with `CCI_scores` populated (a
wide data frame – rows are `"Sender->Receiver"` group pairs, columns are
LR pairs).

`runCCI.default`: a data frame with `"Sender->Receiver"` row names and
one column per significant LR pair containing the interaction score
`0.5 * log2(receiver * sender + 1)`.

## Methods (by class)

- `runCCI(blisa)`: Method for a `blisa` object. If `CCI_scores` are
  already present and `overwrite = FALSE` (the default), the object is
  returned unchanged. Set `overwrite = TRUE` with a `counts_by_group` to
  recompute and replace existing scores. If no scores exist,
  `counts_by_group` must be supplied and scores are computed and
  attached to `x$CCI_scores`.

- `runCCI(default)`: Default method. Performs the raw CCI computation
  and returns only the scores data frame. Typically called internally;
  use `runCCI.blisa` to compute and attach scores to a `blisa` object in
  one step.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")

# CCI is computed automatically when group is supplied to blisa();
# access the scores directly:
head(result$CCI_scores)

# Or compute / recompute scores explicitly:
binned <- hexBinCells(
  as.data.frame(SpatialExperiment::spatialCoords(spe)),
  SummarizedExperiment::assay(spe, "counts"),
  bin_size = 50, group = spe$cell_type
)
result2 <- runCCI(result, counts_by_group = binned$counts_by_group,
                  overwrite = TRUE)
} # }
```
