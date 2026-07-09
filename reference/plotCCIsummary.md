# Sender-by-receiver heatmap of aggregated CCI scores across LR pairs

Generic function. Aggregates CCI scores across all (or the top-ranked)
ligand-receptor pairs and draws a clustered receiver-by-sender heatmap,
one cell per Sender \\\rightarrow\\ Receiver combination.

## Usage

``` r
plotCCIsummary(x, ...)

# S3 method for class 'blisa'
plotCCIsummary(
  x,
  top_lr = NULL,
  sender = NULL,
  receiver = NULL,
  agg_fun = sum,
  main = NULL,
  ...
)

# Default S3 method
plotCCIsummary(
  x,
  top_lr = NULL,
  sender = NULL,
  receiver = NULL,
  agg_fun = sum,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A `blisa` object or a CCI scores data frame (the `CCI_scores` slot of
  a `blisa` object).

- ...:

  Additional arguments passed to the relevant method.

- top_lr:

  Integer or `NULL`. Number of top-ranked LR pairs (by `sig_numbers`) to
  include before aggregating. LR pairs in `CCI_scores` are already
  ordered by rank, so this takes the first `top_lr` columns. `NULL`
  (default) uses all pairs.

- sender:

  Character vector or `NULL`. If provided, only rows where `Sender` is
  in this vector are kept (AND logic with `receiver`). Default `NULL`
  (all senders).

- receiver:

  Character vector or `NULL`. If provided, only rows where `Receiver` is
  in this vector are kept (AND logic with `sender`). Default `NULL` (all
  receivers).

- agg_fun:

  Function used to aggregate scores across LR pairs for each Sender
  \\\rightarrow\\ Receiver combination. Receives a numeric vector with
  `NA`s already removed. Default `sum`.

- main:

  Character or `NULL`. Title drawn above the heatmap. When supplied, the
  heatmap is drawn with this overall title (via
  [`ComplexHeatmap::draw`](https://rdrr.io/pkg/ComplexHeatmap/man/draw-dispatch.html));
  the `Heatmap` object is returned invisibly. Default `NULL` (no title;
  object returned for the caller to print).

## Value

A `Heatmap` object.

## Methods (by class)

- `plotCCIsummary(blisa)`: Method for a `blisa` object. Stops with an
  informative error if `CCI_scores` is `NULL`.

- `plotCCIsummary(default)`: Method for a CCI scores data frame (e.g.
  the `CCI_scores` slot of a `blisa` object).

## See also

[`plotCCILR`](https://chenlaboratory.github.io/blisa/reference/plotCCILR.md)
for a per-LR-pair version of this plot;
[`plotCCI`](https://chenlaboratory.github.io/blisa/reference/plotCCI.md)
for a heatmap with LR pairs as columns.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")
plotCCIsummary(result)
plotCCIsummary(result, top_lr = 10, agg_fun = mean)
} # }
```
