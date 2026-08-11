# Aggregate BLISA hotspots to the pathway level

Rolls up per-ligand-receptor (LR) pair BLISA hotspot results to the
signalling *pathway* level, grouping LR pairs by a pathway annotation
column (e.g. CellChatDB's `pathway_name`). For each pathway it returns
the spots that are hotspots for the pathway, a per-spot pathway p-value,
and the LR pair driving each spot.

## Usage

``` r
blisaPathway(
  x,
  method = c("minP", "simes", "fisher"),
  pathway_col = "pathway_name",
  p_cutoff = 0.05,
  p_adjust = c("none", "BH")
)
```

## Arguments

- x:

  A `blisa` object (or its `LR_results` data frame).

- method:

  Character. Per-spot combination: `"minP"` (default, the original
  behaviour), `"simes"`, or `"fisher"`.

- pathway_col:

  Character. Column of `LR_results` holding the pathway label. Default
  `"pathway_name"` (CellChatDB). Rows with `NA` are dropped.

- p_cutoff:

  Numeric. Significance threshold for `"simes"`/ `"fisher"`. Default
  `0.05`. Ignored for `"minP"`.

- p_adjust:

  Character. Multiple-testing adjustment across spots for
  `"simes"`/`"fisher"`: `"none"` (default) or `"BH"`. Ignored for
  `"minP"`.

## Value

A `blisa` object whose `LR_results` slot is the pathway-level table (one
row per pathway, row names = pathway, ordered by `sig_numbers`
descending), so the standard blisa methods work directly
([`plotHotspots`](https://chenlaboratory.github.io/blisa/reference/plotHotspots.md),
[`plotLRrank`](https://chenlaboratory.github.io/blisa/reference/plotLRrank.md),
`print`). The `bins` and `spatial_weights` slots are carried over from
`x` (when it is a `blisa` object); `CCI_scores` is `NULL`. The
`LR_results` columns are:

- sig_numbers:

  Number of pathway hotspot spots.

- n_LR_pairs:

  Number of LR pairs annotated to the pathway.

- sig_index:

  List column: the hotspot spot indices.

- sig_pval:

  List column: the per-spot pathway p-values.

- LR_pairs:

  List column: all LR pair IDs in the pathway.

- top_LR_pair:

  List column: per spot, the LR pair(s) with the smallest contributing
  p-value.

## Details

The per-spot pathway p-value is computed by one of three `method`s:

`"minP"`

:   

`"simes"`

:   Per-spot Simes combination across *all* LR pairs in the pathway,
    using the full per-spot p-values (`all_pval`) restricted to the
    High-High co-expression direction. Non-High-High/untested
    contributions are set to 1.

`"fisher"`

:   As `"simes"` but Fisher's method.

For `"simes"`/`"fisher"`, a spot is a pathway hotspot when its combined
p-value is \\\le\\ `p_cutoff` (after optional `p_adjust`).

## See also

[`blisa`](https://chenlaboratory.github.io/blisa/reference/blisa.md)

## Examples

``` r
if (FALSE) { # \dontrun{
res <- blisa(spe, platform = "visium")
blisaPathway(res)                          # min-P (original behaviour)
blisaPathway(res, method = "simes", p_adjust = "BH")
} # }
```
