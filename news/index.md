# Changelog

## blisa 1.0.0

- [`blisa()`](https://chenlaboratory.github.io/blisa/reference/blisa.md)
  now uses the `fastLISA` C/OpenMP backend by default (`fast = TRUE`)
  for the bivariate local Moran’s I computation. `fastLISA` moved from
  Suggests to Imports. Results are statistically equivalent to the
  previous
  [`spdep::localmoran_bv`](https://r-spatial.github.io/spdep/reference/localmoran_bv.html)
  backend but not bit-identical (different RNG).
- New `min_total_counts` argument to
  [`blisa()`](https://chenlaboratory.github.io/blisa/reference/blisa.md)
  and
  [`hexBinCells()`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  drops bins whose total counts (summed over all genes) fall below the
  threshold during binning. Default `10`.
- New `verbose` argument (default `FALSE`) to
  [`blisa()`](https://chenlaboratory.github.io/blisa/reference/blisa.md)
  and
  [`hexBinCells()`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  gates progress messages;
  [`blisa()`](https://chenlaboratory.github.io/blisa/reference/blisa.md)
  also shows a progress bar over the ligand-receptor pairs in
  interactive sessions.
- New `lr_pairs` argument to
  [`plotCCI()`](https://chenlaboratory.github.io/blisa/reference/plotCCI.md)
  selects specific ligand-receptor pairs to display, overriding
  `top_lr`.
- New `main` argument to
  [`plotCCI()`](https://chenlaboratory.github.io/blisa/reference/plotCCI.md),
  [`plotCCILR()`](https://chenlaboratory.github.io/blisa/reference/plotCCILR.md),
  and
  [`plotCCIsummary()`](https://chenlaboratory.github.io/blisa/reference/plotCCIsummary.md)
  adds a title above the heatmap.

## blisa 0.2.3

- `LR_results` row names are now informative IDs built from the
  ligand/receptor symbols, with subunits joined by `|` and the two sides
  by `_` (e.g. `TGFB1_TGFBR2|TGFBR1`). The original CellChatDB
  `interaction_name` is retained in a new `interaction_name` column.
- [`plotCCIspatial()`](https://chenlaboratory.github.io/blisa/reference/plotCCIspatial.md)
  now handles multi-subunit ligand/receptor complexes (previously
  errored with `subscript out of bounds`).
- [`plotHotspots()`](https://chenlaboratory.github.io/blisa/reference/plotHotspots.md),
  [`plotCCIspatial()`](https://chenlaboratory.github.io/blisa/reference/plotCCIspatial.md),
  and
  [`plotCCILR()`](https://chenlaboratory.github.io/blisa/reference/plotCCILR.md)
  now match the `ligand`/`receptor` arguments by subunit set,
  independent of order, separator (`,`, `_`, `|`), and whitespace.

## blisa 0.2.2

- Added `fast = TRUE` option to
  [`blisa()`](https://chenlaboratory.github.io/blisa/reference/blisa.md)
  which uses
  [`fastLISA::local_moran_bv`](https://rdrr.io/pkg/fastLISA/man/local_moran_bv.html)
  (C/OpenMP backend) instead of
  [`spdep::localmoran_bv`](https://r-spatial.github.io/spdep/reference/localmoran_bv.html)
  for the bivariate local Moran’s I computation. Default behaviour
  unchanged.
- New `cpu_threads` argument controlling parallelism when `fast = TRUE`.

## blisa 0.2.1

- Update vignettes.

## blisa 0.2.0

CRAN release: 2026-06-02

- Initial CRAN release.
