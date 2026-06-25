filter_genes <- function(gene_column, gene_list) { # filter out the whole interaction if missing any L/R subunit
  # Use parse_units() so subunit splitting matches get_min_expr() exactly
  # (tolerant of ",", "_", and whitespace); the two stages can't diverge.
  vapply(as.character(gene_column), function(s) {
    genes <- parse_units(s)
    if (all(genes %in% gene_list)) {
      paste(genes, collapse = ", ")   # keep: canonical ", " form
    } else {
      NA_character_                   # drop if any subunit is off-panel
    }
  }, character(1), USE.NAMES = FALSE)
}


getLRpairs <- function(gene_panel, LR_df = NULL, species = c("human", "mouse")){
  if (is.null(LR_df)) LR_df <- .load_cellchat_db(species)

  LR_df$ligand.symbol <- filter_genes(LR_df$ligand.symbol, gene_panel)
  LR_df$receptor.symbol <- filter_genes(LR_df$receptor.symbol, gene_panel)

  LR_included <- LR_df[!(is.na(LR_df$ligand.symbol) | is.na(LR_df$receptor.symbol)), ]

  return(LR_included)
}

# NOTE: use ligand/receptor.symbol as gene name, ligand/receptor columns are bilogical entity names which may not be a gene, e.g. RA-ALDH1A3
#' Filter ligand-receptor pairs by expression threshold
#'
#' Retains only LR pairs where at least one bin/spot has counts at or above
#' \code{min_ligand} for every ligand subunit and \code{min_receptor} for every
#' receptor subunit.
#'
#' @param counts Gene-by-bin count matrix (dense or sparse). Row names must be
#'   gene symbols.
#' @param min_ligand Numeric. Minimum count threshold for ligand genes. At least
#'   one bin must meet or exceed this value. Default 10.
#' @param min_receptor Numeric. Minimum count threshold for receptor genes. At
#'   least one bin must meet or exceed this value. Default 10.
#' @param LR_df Data frame of ligand-receptor pairs with columns
#'   \code{ligand.symbol} and \code{receptor.symbol} (comma-separated gene
#'   symbols for multi-subunit complexes). When \code{NULL}, the CellChatDB for
#'   the chosen \code{species} is downloaded automatically.
#' @param species Character. Which CellChatDB to download when \code{LR_df} is
#'   \code{NULL}. One of \code{"human"} (default) or \code{"mouse"}.
#'
#' @return A subset of \code{LR_df} containing only pairs that pass the
#'   expression thresholds for both ligand and receptor.
#' @examples
#' \dontrun{
#' # Supply a small custom LR_df to avoid a network download
#' LR_df <- data.frame(
#'   ligand.symbol   = c("GENE1", "GENE3"),
#'   receptor.symbol = c("GENE2", "GENE4"),
#'   annotation      = c("Secreted Signaling", "ECM-Receptor"),
#'   row.names       = c("LR1", "LR2")
#' )
#' set.seed(1)
#' counts <- matrix(
#'   rpois(4 * 50, lambda = c(20, 1, 5, 20)), nrow = 4, ncol = 50,
#'   dimnames = list(c("GENE1", "GENE2", "GENE3", "GENE4"),
#'                   paste0("bin_", 1:50))
#' )
#' filterLRpairs(counts, min_ligand = 10, min_receptor = 10, LR_df = LR_df)
#' }
#' @export
filterLRpairs <- function(counts, min_ligand = 10, min_receptor = 10,
                          LR_df = NULL, species = c("human", "mouse")) {
  if (is.null(LR_df)) LR_df <- .load_cellchat_db(species)

  filtered_counts_ligand <- counts[Matrix::rowSums(counts >= min_ligand) > 0, ] # at least one spot has more than min_ligand counts
  filtered_counts_receptor <- counts[Matrix::rowSums(counts >= min_receptor) > 0, ]

  LR_filtered_ligand <- getLRpairs(gene_panel = rownames(filtered_counts_ligand), LR_df)
  LR_filtered_receptor <- getLRpairs(gene_panel = rownames(filtered_counts_receptor), LR_df)

  intersect_LR <- intersect(rownames(LR_filtered_ligand), rownames(LR_filtered_receptor))

  LR_filtered <- LR_filtered_ligand[intersect_LR,]

  return(LR_filtered)

}

