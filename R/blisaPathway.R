# Combine a spots-by-LR-pair matrix of p-values into one p-value per spot.
# P is already masked so non-signalling (non-High-High / untested) entries are 1.
.combine_pvals <- function(P, method) {
  k <- ncol(P)
  P <- pmax(P, .Machine$double.xmin)          # guard against log(0)
  switch(method,
    simes = apply(P, 1, function(p) {
      ps <- sort(p)
      min(k * ps / seq_len(k))                # Simes (valid under PRDS)
    }),
    fisher = {
      stat <- -2 * rowSums(log(P))            # Fisher (assumes independence)
      stats::pchisq(stat, df = 2 * k, lower.tail = FALSE)
    }
  )
}

#' Aggregate BLISA hotspots to the pathway level
#'
#' Rolls up per-ligand-receptor (LR) pair BLISA hotspot results to the signalling
#' \emph{pathway} level, grouping LR pairs by a pathway annotation column (e.g.
#' CellChatDB's \code{pathway_name}). For each pathway it returns the spots that
#' are hotspots for the pathway, a per-spot pathway p-value, and the LR pair
#' driving each spot.
#'
#' @details
#' The per-spot pathway p-value is computed by one of three \code{method}s:
#' \describe{
#'   \item{\code{"minP"}}
#'     {The union of each LR pair's significant hotspot spots, with the minimum
#'     p-value per spot. \code{p_cutoff} and \code{p_adjust} are ignored.}
#'   \item{\code{"simes"}}{Per-spot Simes combination across \emph{all} LR pairs
#'     in the pathway, using the full per-spot p-values (\code{all_pval})
#'     restricted to the High-High co-expression direction.
#'     Non-High-High/untested contributions are set to 1.}
#'   \item{\code{"fisher"}}{As \code{"simes"} but Fisher's method.}
#' }
#' For \code{"simes"}/\code{"fisher"}, a spot is a pathway hotspot when its
#' combined p-value is \eqn{\le} \code{p_cutoff} (after optional \code{p_adjust}).
#'
#' @param x A \code{blisa} object (or its \code{LR_results} data frame).
#' @param method Character. Per-spot combination: \code{"minP"} (default,
#'   the original behaviour), \code{"simes"}, or \code{"fisher"}.
#' @param pathway_col Character. Column of \code{LR_results} holding the pathway
#'   label. Default \code{"pathway_name"} (CellChatDB). Rows with \code{NA} are
#'   dropped.
#' @param p_cutoff Numeric. Significance threshold for \code{"simes"}/
#'   \code{"fisher"}. Default \code{0.05}. Ignored for \code{"minP"}.
#' @param p_adjust Character. Multiple-testing adjustment across spots for
#'   \code{"simes"}/\code{"fisher"}: \code{"none"} (default) or \code{"BH"}.
#'   Ignored for \code{"minP"}.
#'
#' @return A \code{blisa} object whose \code{LR_results} slot is the pathway-level
#'   table (one row per pathway, row names = pathway, ordered by
#'   \code{sig_numbers} descending), so the standard blisa methods work directly
#'   (\code{\link{plotHotspots}}, \code{\link{plotLRrank}}, \code{print}). The
#'   \code{bins} and \code{spatial_weights} slots are carried over from \code{x}
#'   (when it is a \code{blisa} object); \code{CCI_scores} is \code{NULL}. The
#'   \code{LR_results} columns are:
#'   \describe{
#'     \item{sig_numbers}{Number of pathway hotspot spots.}
#'     \item{n_LR_pairs}{Number of LR pairs annotated to the pathway.}
#'     \item{sig_index}{List column: the hotspot spot indices.}
#'     \item{sig_pval}{List column: the per-spot pathway p-values.}
#'     \item{LR_pairs}{List column: all LR pair IDs in the pathway.}
#'     \item{top_LR_pair}{List column: per spot, the LR pair(s) with the smallest
#'       contributing p-value.}
#'   }
#'
#' @examples
#' \dontrun{
#' res <- blisa(spe, platform = "visium")
#' blisaPathway(res)                          # min-P (original behaviour)
#' blisaPathway(res, method = "simes", p_adjust = "BH")
#' }
#' @seealso \code{\link{blisa}}
#' @export
blisaPathway <- function(x,
                         method      = c("minP", "simes", "fisher"),
                         pathway_col = "pathway_name",
                         p_cutoff    = 0.05,
                         p_adjust    = c("none", "BH")) {
  method   <- match.arg(method)
  p_adjust <- match.arg(p_adjust)

  LR <- if (inherits(x, "blisa")) x$LR_results else x
  if (!pathway_col %in% colnames(LR))
    stop("Column '", pathway_col, "' not found in LR_results. CellChatDB ",
         "provides 'pathway_name'; set 'pathway_col' if your annotation differs.")
  if (method != "minP" &&
      !all(c("all_pval", "all_quadrant") %in% colnames(LR)))
    stop("'", method, "' needs 'all_pval' and 'all_quadrant' in LR_results; ",
         "ensure blisa() was run to completion.")

  pw   <- LR[[pathway_col]]
  keep <- !is.na(pw)
  LR   <- LR[keep, , drop = FALSE]
  pw   <- pw[keep]
  groups <- split(seq_len(nrow(LR)), pw)

  per_pathway <- lapply(groups, function(rows) {
    lr_ids <- rownames(LR)[rows]

    if (method == "minP") {
      # Union of per-LR significant spots; minimum p-value per spot.
      idx  <- unlist(LR$sig_index[rows])
      pval <- unlist(LR$sig_pval[rows])
      lrv  <- rep(lr_ids, lengths(LR$sig_index[rows]))
      if (length(idx) == 0L)
        return(list(n = 0L, index = integer(0), pval = numeric(0),
                    lr_ids = lr_ids, top = list()))
      min_by_spot <- tapply(pval, idx, min, na.rm = TRUE)
      spots <- as.integer(names(min_by_spot))
      cp    <- as.numeric(min_by_spot)
      top   <- lapply(spots, function(s) {
        here <- idx == s
        lrv[here][pval[here] == min(pval[here], na.rm = TRUE)]
      })
    } else {
      # Proper per-spot combination across ALL LR pairs in the pathway, using
      # the full per-spot p-values restricted to the High-High direction.
      P <- do.call(cbind, LR$all_pval[rows])
      Q <- do.call(cbind, LR$all_quadrant[rows])
      P[is.na(Q) | Q != "High-High"] <- 1
      cp_all <- .combine_pvals(P, method)
      if (p_adjust == "BH") cp_all <- stats::p.adjust(cp_all, "BH")
      spots <- which(cp_all <= p_cutoff)
      cp    <- cp_all[spots]
      top   <- lapply(spots, function(s) lr_ids[which.min(P[s, ])])
    }

    ord <- order(cp)                          # report best spots first
    list(n = length(spots), index = spots[ord], pval = cp[ord],
         lr_ids = lr_ids, top = top[ord])
  })

  # Build a pathway-level "LR_results" table using the standard blisa column
  # names (sig_numbers/sig_index/sig_pval) so the existing blisa methods
  # (plotHotspots, plotLRrank, print) work directly. Rows are pathways.
  LR_out <- do.call(rbind, lapply(names(per_pathway), function(nm) {
    r <- per_pathway[[nm]]
    data.frame(pathway     = nm,
               sig_numbers = r$n,
               n_LR_pairs  = length(r$lr_ids),
               sig_index   = I(list(r$index)),
               sig_pval    = I(list(r$pval)),
               LR_pairs    = I(list(r$lr_ids)),
               top_LR_pair = I(list(r$top)),
               stringsAsFactors = FALSE)
  }))
  LR_out <- as.data.frame(LR_out)
  rownames(LR_out) <- LR_out$pathway
  LR_out$pathway <- NULL
  LR_out <- LR_out[order(-LR_out$sig_numbers), , drop = FALSE]

  # Wrap as a blisa object, carrying the spatial bins/weights when available so
  # plotHotspots() can render pathway hotspots. CCI_scores is not meaningful at
  # the pathway level.
  bins <- if (inherits(x, "blisa")) x$bins            else NULL
  sw   <- if (inherits(x, "blisa")) x$spatial_weights else NULL
  new_blisa(LR_out, bins, sw, CCI_scores = NULL, level = "pathway")
}
