# Run BLISA spatial cell-cell communication analysis

Generic function for running BLISA (Bivariate Local Indicator of Spatial
Association). Dispatches on the class of `x`:

- `blisa.default` accepts a pre-binned gene-by-bin count matrix and a
  matching `bins` polygon object.

- `blisa.SpatialExperiment` accepts a cell-level `SpatialExperiment`
  object and bins cells into hexagonal tiles internally via
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  before running the analysis.

## Usage

``` r
blisa(x, ...)

# Default S3 method
blisa(
  x,
  bins,
  LR_df = NULL,
  bin_size = 50,
  dmax = 250,
  nsim = 999,
  p_cutoff = 0.05,
  min_ligand = 10,
  min_receptor = 10,
  min_cells = 1,
  n_cells_col = NA,
  annotation_col = "annotation",
  default_mode = "diffuse",
  diffuse_category = c("Secreted Signaling", "Non-protein Signaling"),
  species = c("human", "mouse"),
  genes = NULL,
  counts_by_group = NULL,
  fast = TRUE,
  cpu_threads = 4L,
  verbose = FALSE,
  ...
)

# S3 method for class 'SpatialExperiment'
blisa(
  x,
  bin_size = 50,
  LR_df = NULL,
  group = "cell_type",
  genes = NULL,
  min_cells = 1,
  min_total_counts = 10,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A gene-by-bin count matrix (for `blisa.default`) or a cell-level
  `SpatialExperiment` object (for `blisa.SpatialExperiment`).

- ...:

  Additional arguments passed to the relevant method.

- bins:

  An `sf` object of bin polygons. Row order must match the column order
  of `x`.

- LR_df:

  Data frame of ligand-receptor pairs with columns `ligand.symbol` and
  `receptor.symbol`. When `NULL`, CellChatDB for the chosen `species` is
  downloaded automatically.

- bin_size:

  Numeric. Width of each hexagonal bin in coordinate units (e.g.
  microns). Passed to
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  and
  [`computeSpatialWeights`](https://chenlaboratory.github.io/blisa/reference/computeSpatialWeights.md).
  Default `50`.

- dmax:

  Numeric. Maximum distance for diffuse-mode neighbours. Default `250`.

- nsim:

  Integer. Number of permutations for Moran's I significance. Default
  `999`.

- p_cutoff:

  Numeric. P-value threshold for High-High hotspots. Default `0.05`.

- min_ligand:

  Numeric. Minimum ligand count threshold. Default `10`.

- min_receptor:

  Numeric. Minimum receptor count threshold. Default `10`.

- min_cells:

  Integer. Bins with fewer cells are dropped during binning by
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md).
  Default `1`.

- n_cells_col:

  Character or `NA`. Column in `bins` holding per-bin cell counts used
  for `min_cells` filtering. Set to `NA` to skip (default).

- annotation_col:

  Character. Column in `LR_df` specifying interaction category used for
  communication-mode assignment. Default `"annotation"`.

- default_mode:

  Character. CCC mode assigned to LR pairs whose annotation does not
  match `diffuse_category`. Default `"diffuse"`.

- diffuse_category:

  Character vector of annotation categories treated as diffuse
  signalling.

- species:

  Character. Which CellChatDB to download when `LR_df = NULL`. One of
  `"human"` (default) or `"mouse"`.

- genes:

  Character vector of gene names to consider when matching
  ligand-receptor pairs. Defaults to `rownames(x)` (all genes in the
  `SpatialExperiment` object).

- counts_by_group:

  Named list of gene-by-bin count matrices, one per group level (e.g.
  cell type), as returned by
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  when `group` is supplied. When provided,
  [`runCCI`](https://chenlaboratory.github.io/blisa/reference/runCCI.md)
  is called automatically after the BLISA loop and its output is
  included in the result as `CCI_scores`. Default `NULL`.

- fast:

  Logical. When `TRUE` (default), uses
  [`fastLISA::local_moran_bv`](https://rdrr.io/pkg/fastLISA/man/local_moran_bv.html)
  (a fast C/OpenMP backend) for the bivariate local Moran's I
  computation. When `FALSE`, uses the original
  [`spdep::localmoran_bv`](https://r-spatial.github.io/spdep/reference/localmoran_bv.html) +
  [`spdep::hotspot`](https://r-spatial.github.io/spdep/reference/hotspotmap.html)
  pipeline.

- cpu_threads:

  Integer. Number of OpenMP threads used by
  [`fastLISA::local_moran_bv`](https://rdrr.io/pkg/fastLISA/man/local_moran_bv.html).
  Only used when `fast = TRUE`. Ignored on platforms without OpenMP.
  Default `4L`.

- verbose:

  Logical. If `TRUE`, print progress messages. Default `FALSE` (silent).

- group:

  Character. Column name in `colData(x)` to use as the grouping variable
  (e.g. cell type) for per-group bin aggregation and downstream CCI
  analysis via
  [`runCCI`](https://chenlaboratory.github.io/blisa/reference/runCCI.md).
  If the column is not found in `colData(x)`, a message is issued and
  CCI is skipped. Default `"cell_type"`.

- min_total_counts:

  Numeric. Bins whose total counts (summed over all genes) fall below
  this threshold are dropped during binning by
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md).
  Set to `0` to disable. Default `10`.

## Value

A list; see individual method documentation for details.

An object of class `blisa` with four components:

- LR_results:

  Data frame of BLISA results for each LR pair, including `ccc_mode`,
  `sig_numbers`, `sig_index`, `sig_pval`, `all_pval`, `all_lisa`,
  `all_quadrant`, and original columns from `LR_df`. `all_quadrant` is a
  character vector of hotspot quadrant labels (`"High-High"`,
  `"Low-Low"`, etc.) for every bin; non-tested bins are `NA`.

- bins:

  Bin-level `sf` object of hexagonal polygons.

- spatial_weights:

  Spatial weights list from
  [`computeSpatialWeights`](https://chenlaboratory.github.io/blisa/reference/computeSpatialWeights.md).

- CCI_scores:

  `NULL` unless `counts_by_group` is supplied, in which case a wide data
  frame of interaction scores from
  [`runCCI`](https://chenlaboratory.github.io/blisa/reference/runCCI.md):
  rows are `"Sender->Receiver"` group pairs, columns are LR pairs.

## Methods (by class)

- `blisa(default)`: Method for a gene-by-bin count matrix.

- `blisa(SpatialExperiment)`: Method for a cell-level SpatialExperiment
  object. Bins cells into hexagonal tiles via
  [`hexBinCells`](https://chenlaboratory.github.io/blisa/reference/hexBinCells.md)
  then delegates to `blisa.default`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Download and cache the example Xenium breast cancer dataset (~21 MB)
data_url  <- paste0(
  "https://github.com/ChenLaboratory/example_data/releases/",
  "download/v1.0.0/spe_xenium_bc_s1rep1.rds"
)
cache_dir <- tools::R_user_dir("blisa", "cache")
data_file <- file.path(cache_dir, "spe_xenium_bc_s1rep1.rds")
if (!file.exists(data_file)) {
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  download.file(data_url, data_file, mode = "wb")
}
spe <- readRDS(data_file)

# blisa.SpatialExperiment method: bins cells, runs LISA, scores CCI
result <- blisa(spe, bin_size = 50, group = "cell_type")
result
} # }
```
