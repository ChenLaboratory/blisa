# Suppress R CMD check notes for ggplot2 aes() variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # ggplot2 aes (plotCCIspatial)
    "cell_pair_plot",
    # ggplot2 aes (plotHotspots)
    "fill_col", "val", "y",
    # ggplot2 aes (plotLRrank)
    "LR_pair", "sig_numbers", "annotation"
  ))
}

LRI_spatial_colors <- c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")

col.pDark <- c(
  "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
)

col.pLight <- c(
  "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#C7C7C7", "#DBDB8D", "#9EDAE5"
)

col.pMedium <- c(
  "#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
  "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA"
)

cols <- c(col.pDark, col.pLight, col.pMedium)

# Session-level cache so CellChatDB is only downloaded once per R session
.blisa_cache <- new.env(parent = emptyenv())

.load_cellchat_db <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  key <- paste0("CellChatDB.", species)
  if (!is.null(.blisa_cache[[key]])) return(.blisa_cache[[key]])
  url <- sprintf(
    "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.%s.rda",
    species
  )
  message("Downloading CellChatDB.", species, " from GitHub (once per session)...")
  tmp_file <- tempfile(fileext = ".rda")
  on.exit(unlink(tmp_file), add = TRUE)
  tryCatch({
    download.file(url, tmp_file, mode = "wb", quiet = TRUE)
    tmp <- new.env(parent = emptyenv())
    load(tmp_file, envir = tmp)
    .blisa_cache[[key]] <- tmp[[key]]$interaction
  }, error = function(e) stop(
    "Failed to download CellChatDB.", species, ". ",
    "Check your internet connection or supply LR_df explicitly."
  ))
  .blisa_cache[[key]]
}

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
