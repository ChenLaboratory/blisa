filter_genes <- function(gene_column, gene_list) { # function to filter out the whole interaction if missing any L/R subunits
  sapply(strsplit(as.character(gene_column), ", "), function(genes) {
    filtered_genes <- genes[genes %in% gene_list]
    if (length(filtered_genes) != length(genes)) {
      return(NA)
    } else {
      return(paste(filtered_genes, collapse = ", "))
    }
  })
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

# LR_pairs_filtered <- filterLRpairs(counts = hex_gene_counts, min_ligand = 10, min_receptor = 10, LR_df = CellChatDB.human$interaction)
