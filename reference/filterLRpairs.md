# Filter ligand-receptor pairs by expression threshold

Retains only LR pairs where at least one bin/spot has counts at or above
`min_ligand` for every ligand subunit and `min_receptor` for every
receptor subunit.

## Usage

``` r
filterLRpairs(
  counts,
  min_ligand = 10,
  min_receptor = 10,
  LR_df = NULL,
  species = c("human", "mouse")
)
```

## Arguments

- counts:

  Gene-by-bin count matrix (dense or sparse). Row names must be gene
  symbols.

- min_ligand:

  Numeric. Minimum count threshold for ligand genes. At least one bin
  must meet or exceed this value. Default 10.

- min_receptor:

  Numeric. Minimum count threshold for receptor genes. At least one bin
  must meet or exceed this value. Default 10.

- LR_df:

  Data frame of ligand-receptor pairs with columns `ligand.symbol` and
  `receptor.symbol` (comma-separated gene symbols for multi-subunit
  complexes). When `NULL`, the CellChatDB for the chosen `species` is
  downloaded automatically.

- species:

  Character. Which CellChatDB to download when `LR_df` is `NULL`. One of
  `"human"` (default) or `"mouse"`.

## Value

A subset of `LR_df` containing only pairs that pass the expression
thresholds for both ligand and receptor.

## Examples

``` r
if (FALSE) { # \dontrun{
# Supply a small custom LR_df to avoid a network download
LR_df <- data.frame(
  ligand.symbol   = c("GENE1", "GENE3"),
  receptor.symbol = c("GENE2", "GENE4"),
  annotation      = c("Secreted Signaling", "ECM-Receptor"),
  row.names       = c("LR1", "LR2")
)
set.seed(1)
counts <- matrix(
  rpois(4 * 50, lambda = c(20, 1, 5, 20)), nrow = 4, ncol = 50,
  dimnames = list(c("GENE1", "GENE2", "GENE3", "GENE4"),
                  paste0("bin_", 1:50))
)
filterLRpairs(counts, min_ligand = 10, min_receptor = 10, LR_df = LR_df)
} # }
```
