# Changelog

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
