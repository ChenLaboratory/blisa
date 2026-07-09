# Heatmap of CCI scores across all ligand-receptor pairs

Generic function. Draws a clustered heatmap (via `ComplexHeatmap`) with
rows as Sender \\\rightarrow\\ Receiver cell-type pairs and columns as
LR pairs. Row annotations colour-code the sender and receiver cell
types.

## Usage

``` r
plotCCI(x, ...)

# S3 method for class 'blisa'
plotCCI(
  x,
  top_lr = 20,
  top_pairs = 30,
  lr_pairs = NULL,
  sender = NULL,
  receiver = NULL,
  colors = NULL,
  colours = NULL,
  main = NULL,
  ...
)

# Default S3 method
plotCCI(
  x,
  top_lr = 20,
  top_pairs = 30,
  lr_pairs = NULL,
  sender = NULL,
  receiver = NULL,
  colors = NULL,
  colours = NULL,
  main = NULL,
  ...
)
```

## Arguments

- x:

  A `blisa` object or a CCI scores data frame (the `CCI_scores` slot of
  a `blisa` object). The data frame must contain columns `Sender`,
  `Receiver`, and one column per LR pair.

- ...:

  Additional arguments passed to the relevant method.

- top_lr:

  Integer or `NULL`. Number of top-ranked LR pairs (by `sig_numbers`) to
  display as columns. LR pairs in `CCI_scores` are already ordered by
  rank, so this simply takes the first `top_lr` columns. Ignored when
  `lr_pairs` is supplied. Default `20`.

- top_pairs:

  Integer or `NULL`. Number of top sender-receiver pairs to display as
  rows, ranked by their maximum interaction score across the displayed
  LR pairs (after `top_lr` / `lr_pairs` is applied). When `NULL` all
  rows are shown. Default `30`.

- lr_pairs:

  Character vector or `NULL`. When supplied, only these LR pair column
  names are shown, in the order given, overriding `top_lr`. Names not
  found in `CCI_scores` are dropped with a warning. Default `NULL` (use
  `top_lr`).

- sender:

  Character vector or `NULL`. If provided, only rows where `Sender` is
  in this vector are kept. Applied independently of `receiver` (AND
  logic when both are supplied). Default `NULL` (all senders).

- receiver:

  Character vector or `NULL`. If provided, only rows where `Receiver` is
  in this vector are kept. Applied independently of `sender` (AND logic
  when both are supplied). Default `NULL` (all receivers).

- colors:

  Named character vector mapping cell-type names to colours, used for
  the sender/receiver row annotations. When `NULL` (default), colours
  are assigned automatically from the package palette. The British
  spelling `colours` is accepted as an alias.

- colours:

  Alias for `colors` (British spelling). Ignored when `colors` is
  supplied.

- main:

  Character or `NULL`. Title drawn above the heatmap (mapped to the
  `column_title` of
  [`ComplexHeatmap::Heatmap`](https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html)).
  Default `NULL` (no title).

## Value

Invisibly returns the `Heatmap` object.

## Methods (by class)

- `plotCCI(blisa)`: Method for a `blisa` object. Extracts `CCI_scores`
  and delegates to `plotCCI.default`. Stops with an informative error if
  `CCI_scores` is `NULL`.

- `plotCCI(default)`: Method for a CCI scores data frame (e.g. the
  `CCI_scores` slot of a `blisa` object).

## See also

[`plotCCILR`](https://chenlaboratory.github.io/blisa/reference/plotCCILR.md)
for a sender-by-receiver heatmap of a single LR pair;
[`plotCCIsummary`](https://chenlaboratory.github.io/blisa/reference/plotCCIsummary.md)
for an aggregated sender-by-receiver heatmap across all LR pairs.

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")
plotCCI(result, top_lr = 20, top_pairs = 30)
} # }
```
