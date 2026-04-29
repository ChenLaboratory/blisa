# =============================================================================
# Expression utilities
# =============================================================================

# Split "A, B" or "A_B" into a character vector of gene symbols
parse_units <- function(s) {
  s <- gsub("\\s+", "", as.character(s))
  unlist(strsplit(s, "[,_]"))
}

# Per-bin minimum expression across all subunits of a (possibly multi-unit) gene
# hex_gene_counts: genes x bins matrix
get_min_expr <- function(gene_str, hex_gene_counts) {
  genes <- parse_units(gene_str)
  genes <- genes[genes %in% rownames(hex_gene_counts)]
  n     <- ncol(hex_gene_counts)
  if (length(genes) == 0) return(rep(0, n))  # gene(s) absent -> 0
  if (length(genes) == 1) return(as.numeric(hex_gene_counts[genes, ]))
  # per-bin minimum across subunits
  mat <- hex_gene_counts[genes, , drop = FALSE]
  apply(mat, 2, min)
}

# Assign ccc_mode ("diffuse" / "nearby") based on annotation column
LR_df_add_mode <- function(LR_df,
                           col              = "annotation",
                           default_mode     = "diffuse",
                           diffuse_category = c("Secreted Signaling",
                                                "Non-protein Signaling")) {
  # Case 1: annotation column missing → all diffuse
  if (!col %in% colnames(LR_df)) {
    LR_df$ccc_mode <- default_mode
    message(col, " column missing — setting ccc_mode='", default_mode, "' for all.")
    return(LR_df)
  }
  
  # Case 2: annotation exists
  LR_df$ccc_mode <-
    ifelse(
      LR_df[[col]] %in% diffuse_category,
      "diffuse",
      "nearby"
    )
  
  message("ccc_mode: 'diffuse' for [",
          paste(diffuse_category, collapse = ", "),
          "]; 'nearby' for others.")
  LR_df
}


#' Run BLISA on Bin-Level Spatial Data with Isolate Removal
#'
#' Performs BLISA (Bivariate Local Indicator of Spatial Association) analysis
#' on bin-level spatial transcriptomics data to identify spatially enriched
#' ligand-receptor interactions. Spatially isolated bins and bins failing
#' minimum cell count filtering are excluded from Moran's I calculation and
#' assigned neutral statistics in the final output (p = 1, I = 0).
#'
#' @param counts_matrix Numeric gene-by-bin expression matrix. Rows correspond
#'   to genes and columns correspond to spatial bins.
#' @param bin_sf An \code{sf} object containing bin-level spatial geometries
#'   and metadata. Must contain one row per bin.
#' @param LR_df Data frame of ligand-receptor pairs. Must contain ligand and
#'   receptor gene symbol columns compatible with \code{filterLRpairs()}.
#' @param hex_size Numeric. Approximate bin/hexagon spacing used to define
#'   queen adjacency for nearby signaling mode.
#' @param dmax Numeric. Maximum diffusion distance (in spatial coordinate units)
#'   used for diffuse signaling mode.
#' @param nsim Integer. Number of permutations for local bivariate Moran's I
#'   significance estimation.
#' @param p_cutoff Numeric. P-value threshold for calling significant
#'   High-High ligand-receptor hotspots.
#' @param min_ligand Numeric. Minimum ligand expression threshold required for
#'   a ligand to be retained in ligand-receptor filtering.
#' @param min_receptor Numeric. Minimum receptor expression threshold required
#'   for a receptor to be retained in ligand-receptor filtering.
#' @param min_cells_per_bin Integer. Minimum number of cells required for a bin
#'   to be included in Moran's I calculation. Bins below this threshold are
#'   excluded and assigned neutral statistics. Ignored when \code{n_cells_col = NA}.
#' @param n_cells_col Character or \code{NA}. Column name in \code{bin_sf}
#'   holding per-bin cell counts used for \code{min_cells_per_bin} filtering.
#'   Set to \code{NA} to skip cell-count filtering (default).
#' @param col Character. Column name in \code{LR_df} specifying ligand-receptor
#'   annotation/category used for communication mode assignment.
#' @param default_mode Character. Default communication mode assigned to ligand-
#'   receptor pairs not matching \code{diffuse_category}. Typically one of
#'   \code{"nearby"} or \code{"diffuse"}.
#' @param diffuse_category Character vector of ligand-receptor annotation
#'   categories to be treated as diffuse signaling.
#'
#' @returns A list containing:
#' \describe{
#'   \item{LR_out}{Data frame summarizing BLISA results for each ligand-receptor
#'   pair, including significant hotspot counts, hotspot indices, p-values,
#'   and local Moran's I statistics.}
#'   \item{bin_sf}{Input bin-level \code{sf} object.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- runBLISA.default.isoaltes.removed(
#'   counts_matrix = counts_mat,
#'   bin_sf = bin_sf,
#'   LR_df = LR_pairs
#' )
#' }
runBLISA.default.isolates.removed <- function(
    counts_matrix,
    bin_sf,
    LR_df,
    hex_size = 50,
    dmax = 250,
    nsim = 999,
    p_cutoff = 0.05,
    min_ligand = 10, min_receptor = 10,
    min_cells_per_bin = 1,
    n_cells_col = NA,
    col = "annotation",
    default_mode = "diffuse",
    diffuse_category = c("Secreted Signaling", "Non-protein Signaling")
) {

  sw <- computeSpatialWeights(bin_sf, hex_size, dmax, min_cells_per_bin, n_cells_col)
  queen_wt <- sw$queen_wt
  dist_wt <- sw$dist_wt
  keep_idx_queen <- sw$keep_idx_queen
  keep_idx_dist <- sw$keep_idx_dist

  ## ---------------------------
  ## Filter LR pairs
  ## ---------------------------
  LR_df_filtered <- filterLRpairs(
    counts = counts_matrix,
    min_ligand,
    min_receptor,
    LR_df
  )

  ## ---------------------------
  ## Local bivariate Moran
  ## ---------------------------
  LR_out <- LR_df_add_mode(
    LR_df_filtered,
    col,
    default_mode,
    diffuse_category
  )

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

    # Perform in subset with no isolates, then restore full-length vectors with neutral values for excluded bins.
    if (mode == "nearby") {
      wt <- queen_wt
      keep_idx <- keep_idx_queen
    } else {
      wt <- dist_wt
      keep_idx <- keep_idx_dist
    }

    message("ccc mode is ", mode)

    x_full <- get_min_expr(receptor, counts_matrix)
    y_full <- get_min_expr(ligand, counts_matrix)

    x <- x_full[keep_idx]
    y <- y_full[keep_idx]

    res_bv <- spdep::localmoran_bv(x, y, wt, nsim = nsim)

    ## Restore full-length vectors
    full_pval <- rep(1, length(x_full))
    full_lisa <- rep(0, length(x_full))

    full_pval[keep_idx] <- res_bv[, "Pr(folded) Sim"]
    full_lisa[keep_idx] <- res_bv[, "Ibvi"]

    hs <- spdep::hotspot(
      res_bv,
      Prname = "Pr(folded) Sim",
      cutoff = p_cutoff,
      quadrant.type = "pysal",
      p.adjust = "none"
    )

    hs_idx_hh <- !is.na(hs) & hs == "High-High"

    HH_idx <- keep_idx[which(hs_idx_hh)]
    HH_pval <- full_pval[HH_idx]

    LR_out$sig_numbers[i] <- length(HH_idx)
    LR_out$sig_index[[i]] <- HH_idx
    LR_out$sig_pval[[i]]  <- HH_pval
    LR_out$all_pval[[i]]  <- full_pval
    LR_out$all_lisa[[i]]  <- full_lisa
    LR_out$ccc_mode[i]    <- mode
  }

  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]

  front_cols <- c("ccc_mode", "sig_numbers", "sig_index", "sig_pval")

  LR_out <- LR_out[, c(front_cols, setdiff(colnames(LR_out), front_cols))]

  list(LR_out = LR_out, bin_sf = bin_sf, sw = sw)
}
