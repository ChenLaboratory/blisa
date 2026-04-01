parse_units <- function(s) {
  s <- as.character(s)
  s <- gsub("\\s+", "", s)
  # allow "A,B" or "A_B" or "A, B"
  unlist(strsplit(s, "[,_]"))
}

# return numeric vector length n_hex
get_min_expr <- function(gene_str, hex_gene_counts) {
  genes <- parse_units(gene_str)
  genes <- genes[genes %in% rownames(hex_gene_counts)]
  if (length(genes) == 0) {
    return(rep(0, n_hex))  # gene(s) absent -> 0
  }
  if (length(genes) == 1) {
    return(as.numeric(hex_gene_counts[genes, ]))
  }
  # per-bin minimum across subunits
  mat <- hex_gene_counts[genes, , drop = FALSE]
  apply(mat, 2, min)
}

LR_df_add_mode <- function(LR_df, col = "annotation", default_mode = "diffuse",
                           diffuse_category = c("Secreted Signaling", "Non-protein Signaling")) {
  # Case 1: annotation column missing → all diffuse
  if (!col %in% colnames(LR_df)) {
    LR_df$ccc_mode <- default_mode
    message(paste0(col, " column is missing. Setting ccc_mode as '", default_mode, "' for all."))

    return(LR_df)
  }

  # Case 2: annotation exists
  LR_df$ccc_mode <-
    ifelse(
      LR_df[[col]] %in% diffuse_category,
      "diffuse",
      "nearby"
    )
  message(paste0("ccc_mode is 'diffuse' for category: ", paste0(diffuse_category, collapse = ", "), "; 'nearby' for others."))

  LR_df
}


hex_binning_cells <- function(
    coords_df = NULL,
    counts_matrix = NULL,
    spe = NULL,
    hex_size = 50) {

  ## ---------------------------
  ## 0. Input handling
  ## ---------------------------
  if (!is.null(spe)) {
    coords_df <- as.data.frame(spe@int_colData$spatialCoords)
    counts_matrix <- spe@assays@data$counts
  } else {
    stopifnot(!is.null(coords_df), !is.null(counts_matrix))
    stopifnot(inherits(counts_matrix, c("matrix", "dgCMatrix")))
  }

  coords_df$cell_id <- rownames(coords_df)

  required_cols <- c("cell_id", "x_centroid", "y_centroid")
  stopifnot(all(required_cols %in% colnames(coords_df)))

  stopifnot(all(colnames(counts_matrix) %in% coords_df$cell_id))

  ## ---------------------------
  ## 1. sf cells
  ## ---------------------------
  cell_sf <- sf::st_as_sf(
    coords_df,
    coords = c("x_centroid", "y_centroid"),
    crs = NA
  )

  ## ---------------------------
  ## 2. Hex binning
  ## ---------------------------
  # create bins
  hex_geom <- sf::st_make_grid( # only polygon info
    cell_sf,
    cellsize = hex_size,
    what = "polygons",
    square = FALSE
  )

  hex_sf <- sf::st_sf( # sf object
    hex_id = seq_along(hex_geom),
    geometry = hex_geom
  )
  dim(hex_sf)

  # map cells to bins
  cell_hex_sf <- sf::st_join(cell_sf, hex_sf, join = sf::st_intersects)
  cell_hex_df <- sf::st_drop_geometry(cell_hex_sf)
  dim(cell_hex_df)
  length(unique(cell_hex_df$hex_id))

  ## ---------------------------
  ## 3. Aggregate counts
  ## ---------------------------
  cell_to_hex <- cell_hex_df$hex_id
  names(cell_to_hex) <- cell_hex_df$cell_id
  cell_to_hex <- cell_to_hex[colnames(counts_matrix)]

  n_hex <- nrow(hex_sf) # hex_sf contains ALL hexes including empty bins
  hex_factor <- factor(cell_to_hex, levels = seq_len(n_hex))

  H <- Matrix::sparse.model.matrix(~ hex_factor - 1)
  hex_gene_counts <- counts_matrix %*% H
  dim(hex_gene_counts)

  hex_sf$n_cells <- as.numeric(Matrix::colSums(H))
  dim(hex_sf)

  list(
    hex_sf = hex_sf,
    hex_gene_counts = hex_gene_counts
  )

}

# add binning step
# runBLISA.default <- function( # bin-level
#     counts_matrix,
#     sf,
#     LR_df, # use index, remove other info
#     hex_size = 50,
#     dmax = 250,
#     nsim = 999,
#     p_cutoff = 0.05,
#     min_ligand = 10, min_receptor = 10,
#     col = "annotation", default_mode = "diffuse",
#     diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
# )
#
# runBLISA.spe <- function(
#     spe, # default cell -> bining inside
#     sf = NULL, # gridding info
#     LR_df = NULL, # default cellchatDB (download to local/ git online read)
#     hex_size = 50,
#     dmax = 250,
#     nsim = 999,
#     p_cutoff = 0.05,
#     min_ligand = 10, min_receptor = 10,
#     col = "annotation", default_mode = "diffuse",
#     diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
# )

#' Run BLISA analysis
#'
#' @param counts_matrix Gene expression matrix
#' @param coords_df spatial coordimates dataframe
#' @param LR_df Ligand-receptor dataframe
#' @export
runBLISA.old <- function(
    coords_df,
    counts_matrix,
    LR_df,
    hex_size = 50,
    dmax = 250,
    nsim = 999,
    p_cutoff = 0.05,
    min_ligand = 10, min_receptor = 10,
    col = "annotation", default_mode = "diffuse",
    diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
) {

  ## ---------------------------
  ## Hex binning
  ## ---------------------------
  binning_res <- hex_binning_cells(coords_df, counts_matrix, spe, hex_size)
  hex_gene_counts <- binning_res$hex_gene_counts
  hex_sf <- binning_res$hex_sf

  ## ---------------------------
  ## Spatial weights
  ## ---------------------------
  centroids <- sf::st_centroid(hex_sf)
  coords <- sf::st_coordinates(centroids)

  # exponential distance weight
  dist_nb <- spdep::dnearneigh(coords, 0, dmax)

  weight_at_dmax <- 0.01 # exp dist weight=0.01 for dist=dmax
  dist_wt <- spdep::nb2listwdist(
    dist_nb,
    hex_sf,
    type = "exp",
    style = "W",
    alpha = -log(weight_at_dmax) / dmax
  )

  # queen weight order = 1
  queen_nb <- spdep::dnearneigh(centroids, 0, 1.2*hex_size)
  queen_wt <- spdep::nb2listwdist(queen_nb, centroids, type="idw", style="W", zero.policy = TRUE)

  ## ---------------------------
  ## Filter LR pairs
  ## ---------------------------
  # only use LR pairs with at least n counts in at least one bin
  LR_df_filtered <- filterLRpairs(counts = hex_gene_counts,
                                  min_ligand, min_receptor,
                                  LR_df)

  ## ---------------------------
  ## Local bivariate Moran for all LR pairs
  ## ---------------------------
  LR_out <- LR_df_add_mode(LR_df_filtered, col, default_mode, diffuse_category)

  LR_out$sig_numbers <- integer(nrow(LR_out))
  LR_out$sig_index   <- vector("list", nrow(LR_out))
  LR_out$sig_pval    <- vector("list", nrow(LR_out))
  LR_out$all_pval    <- vector("list", nrow(LR_out))
  LR_out$all_lisa    <- vector("list", nrow(LR_out))

  for (i in seq_len(nrow(LR_out))) {
    message(rownames(LR_out)[i])

    ligand   <- LR_out$ligand.symbol[i]
    receptor <- LR_out$receptor.symbol[i]

    mode <- LR_out$ccc_mode[i]
    wt <- if (mode == "nearby") queen_wt else dist_wt
    message("ccc mode is ", mode)

    # bivariate vectors per hex (min for multi-unit)
    x <- get_min_expr(receptor, hex_gene_counts)  # receptor
    y <- get_min_expr(ligand, hex_gene_counts)    # ligand

    # bivariate local moran
    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    hs <- spdep::hotspot(
      res_bv,
      Prname = "Pr(folded) Sim",
      cutoff = p_cutoff,
      quadrant.type = "pysal",
      p.adjust = "none"
    )

    idx <- (hs == "High-High")
    idx[is.na(idx)] <- FALSE

    HH_idx  <- which(idx)
    HH_pval <- res_bv[idx, "Pr(folded) Sim"]

    # write back into LR_out
    LR_out$sig_numbers[i] <- length(HH_idx)
    LR_out$sig_index[[i]] <- HH_idx
    LR_out$sig_pval[[i]]  <- HH_pval
    LR_out$all_pval[[i]]  <- res_bv[, "Pr(folded) Sim"]
    LR_out$all_lisa[[i]]  <- res_bv[, "Ibvi"]
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]

  ## ---------------------------
  ## Return both results + binned data
  ## ---------------------------
  front_cols <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")

  LR_out <- LR_out[, c(front_cols, setdiff(colnames(LR_out), front_cols))]

  list(
    LR_out = LR_out,
    hex_sf = hex_sf,
    hex_gene_counts = hex_gene_counts
  )

}



