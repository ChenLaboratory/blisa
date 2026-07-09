# Dot plot ranking LR pairs by number of significant hotspot bins

Generic function for ranking LR pairs. Dispatches on the class of `x`:

- `plotLRrank.blisa` accepts a `blisa` object and uses its `LR_results`
  slot directly.

- `plotLRrank.data.frame` accepts the `LR_results` data frame directly.

## Usage

``` r
plotLRrank(x, ...)

# S3 method for class 'blisa'
plotLRrank(x, top = 30, pt_size = 4, flip = FALSE, ...)

# S3 method for class 'data.frame'
plotLRrank(x, top = 30, pt_size = 4, flip = FALSE, ...)
```

## Arguments

- x:

  A `blisa` object or a data frame of LR results. The data frame must
  contain columns `sig_numbers` and `annotation`.

- ...:

  Additional arguments passed to the relevant method.

- top:

  Integer or `NULL`. Number of top LR pairs (by `sig_numbers`) to
  display. Default `30`.

- pt_size:

  Numeric. Point size passed to `geom_point`. Default 4.

- flip:

  Logical. When `TRUE`, LR pairs are placed on the x-axis and the
  hotspot count on the y-axis (vertical orientation). Default `FALSE`
  (LR pairs on y-axis, horizontal orientation).

## Value

A `ggplot` object.

## Methods (by class)

- `plotLRrank(blisa)`: Method for a `blisa` object. Extracts
  `LR_results` and delegates to `plotLRrank.data.frame`.

- `plotLRrank(data.frame)`: Method for a data frame of LR results (e.g.
  the `LR_results` slot of a `blisa` object).

## Examples

``` r
if (FALSE) { # \dontrun{
# Continuing from the blisa() example:
# result <- blisa(spe, bin_size = 50, group = "cell_type")
plotLRrank(result, top = 30)
plotLRrank(result, top = 20, flip = TRUE)
} # }
```
