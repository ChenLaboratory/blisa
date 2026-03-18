# NOTE: use ligand/receptor.symbol as gene name, ligand/receptor columns are bilogical entity names which may not be a gene, e.g. RA-ALDH1A3

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


getLRpairs <- function(gene_panel, LR_df = CellChatDB.human$interaction){

  LR_df$ligand.symbol <- filter_genes(LR_df$ligand.symbol, gene_panel)
  LR_df$receptor.symbol <- filter_genes(LR_df$receptor.symbol, gene_panel)

  LR_included <- LR_df[!(is.na(LR_df$ligand.symbol) | is.na(LR_df$receptor.symbol)), ]

  return(LR_included)
}

filterLRpairs <- function(counts, min_ligand = 10, min_receptor = 10, LR_df = CellChatDB.human$interaction) {
  #counts = spe@assays@data$counts

  filtered_counts_ligand <- counts[Matrix::rowSums(counts >= min_ligand) > 0, ] # at least one spot has more than min_ligand counts
  filtered_counts_receptor <- counts[Matrix::rowSums(counts >= min_receptor) > 0, ]

  LR_filtered_ligand <- getLRpairs(gene_panel = rownames(filtered_counts_ligand), LR_df)
  LR_filtered_receptor <- getLRpairs(gene_panel = rownames(filtered_counts_receptor), LR_df)

  intersect_LR <- intersect(rownames(LR_filtered_ligand), rownames(LR_filtered_receptor))

  LR_filtered <- LR_filtered_ligand[intersect_LR,]

  return(LR_filtered)

}

# LR_pairs_filtered <- filterLRpairs(counts = hex_gene_counts, min_ligand = 10, min_receptor = 10, LR_df = CellChatDB.human$interaction)
