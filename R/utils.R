LRI_spatial_colors <- c("#FFFFCC", "#FFD700", "#FF7F00", "#D7301F")

cols <- c(
  # Highly distinct / saturated first
  "#E69F00", # orange
  "#56B4E9", # light blue
  "#009E73", # green
  "#D55E00", # vermillion
  "#CC79A7", # magenta
  "#117a77", # blue
  "#F0E442", # yellow
  "#0b2b5e", # strong blue
  "#33A02C", # strong green
  "#E31A1C", # red
  "#6A3D9A", # purple
  "#B15928", # brown

  # Softer but still distinct
  "#FB8072", # coral
  "#80B1D3", # sky blue
  "#FDB462", # orange pastel
  "#B3DE69", # lime pastel
  "#CAB2D6", # lavender
  "#FB9A99", # salmon
  "#a18e6a", # gold pastel
  "#FF7F00", # bright orange

  # Pastel / similar shades pushed later
  "#8DD3C7", # teal pastel
  "#FFFFB3", # pale yellow
  "#BEBADA", # pale lavender
  "#FCCDE5", # pink pastel
  "#BC80BD", # violet pastel
  "#CCEBC5", # mint
  "#FFED6F", # soft yellow
  "#A6CEE3", # pale blue
  "#B2DF8A", # pale green
  "#FFB3BA"  # baby pink
)

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
