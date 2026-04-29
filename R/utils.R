parse_units <- function(s) {
  s <- as.character(s)
  s <- gsub("\\s+", "", s)
  unlist(strsplit(s, "[,_]"))
}

get_min_expr <- function(gene_str, hex_gene_counts) {
  genes <- parse_units(gene_str)
  genes <- genes[genes %in% rownames(hex_gene_counts)]
  n <- ncol(hex_gene_counts)
  if (length(genes) == 0) return(rep(0, n))
  if (length(genes) == 1) return(as.numeric(hex_gene_counts[genes, ]))
  mat <- hex_gene_counts[genes, , drop = FALSE]
  apply(mat, 2, min)
}

LR_df_add_mode <- function(LR_df, col = "annotation", default_mode = "diffuse",
                           diffuse_category = c("Secreted Signaling", "Non-protein Signaling")) {
  if (!col %in% colnames(LR_df)) {
    LR_df$ccc_mode <- default_mode
    message(col, " column is missing. Setting ccc_mode as '", default_mode, "' for all.")
    return(LR_df)
  }
  LR_df$ccc_mode <- ifelse(LR_df[[col]] %in% diffuse_category, "diffuse", "nearby")
  message("ccc_mode is 'diffuse' for category: ", paste(diffuse_category, collapse = ", "), "; 'nearby' for others.")
  LR_df
}
