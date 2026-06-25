# blisa 0.2.3

* plotCCIspatial() now handles multi-subunit ligand/receptor complexes 
  (previously errored with subscript out of bounds).
* New row name formatting for multi-subunit ligand/receptor.

# blisa 0.2.2

* Added `fast = TRUE` option to `blisa()` which uses `fastLISA::local_moran_bv`
  (C/OpenMP backend) instead of `spdep::localmoran_bv` for the bivariate local
  Moran's I computation. Default behaviour unchanged.
* New `cpu_threads` argument controlling parallelism when `fast = TRUE`.

# blisa 0.2.1

* Update vignettes.

# blisa 0.2.0

* Initial CRAN release.
