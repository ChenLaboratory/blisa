
#' Score cell-cell interactions from BLISA hotspots
#'
#' Generic function for scoring cell-cell interactions. Dispatches on the
#' class of \code{x}:
#' \itemize{
#'   \item \code{runCCI.blisa} accepts a \code{blisa} object. If
#'     \code{CCI_scores} are already present and \code{overwrite = FALSE}
#'     (the default), the object is returned unchanged. Set
#'     \code{overwrite = TRUE} with a \code{counts_by_group} to recompute and
#'     replace existing scores. If no scores exist, \code{counts_by_group}
#'     must be supplied and scores are computed and attached.
#'   \item \code{runCCI.default} performs the raw computation given a
#'     \code{blisa} object and a \code{counts_by_group} list, returning only
#'     the scores data frame. Used internally by \code{runCCI.blisa} and
#'     \code{\link{blisa.default}}.
#' }
#'
#' @param x A \code{blisa} object.
#' @param ... Additional arguments passed to the relevant method.
#'
#' @return See individual method documentation.
#' @examples
#' \dontrun{
#' # Continuing from the blisa() example:
#' # result <- blisa(spe, bin_size = 50, group = "cell_type")
#'
#' # CCI is computed automatically when group is supplied to blisa();
#' # access the scores directly:
#' head(result$CCI_scores)
#'
#' # Or compute / recompute scores explicitly:
#' binned <- hexBinCells(
#'   as.data.frame(SpatialExperiment::spatialCoords(spe)),
#'   SummarizedExperiment::assay(spe, "counts"),
#'   bin_size = 50, group = spe$cell_type
#' )
#' result2 <- runCCI(result, counts_by_group = binned$counts_by_group,
#'                   overwrite = TRUE)
#' }
#' @export
runCCI <- function(x, ...) UseMethod("runCCI")


#' @describeIn runCCI Method for a \code{blisa} object. If \code{CCI_scores}
#'   are already present and \code{overwrite = FALSE} (the default), the object
#'   is returned unchanged. Set \code{overwrite = TRUE} with a
#'   \code{counts_by_group} to recompute and replace existing scores. If no
#'   scores exist, \code{counts_by_group} must be supplied and scores are
#'   computed and attached to \code{x$CCI_scores}.
#'
#' @param counts_by_group Named list of gene-by-bin sparse count matrices, one
#'   per group level (e.g. cell type). Typically the \code{counts_by_group}
#'   element of the list returned by \code{\link{hexBinCells}} when
#'   \code{group} is supplied. Names must match the group levels. Required
#'   when \code{x$CCI_scores} is \code{NULL} or when \code{overwrite = TRUE}.
#' @param overwrite Logical. If \code{FALSE} (default) and \code{x$CCI_scores}
#'   is already populated, the object is returned unchanged. If \code{TRUE} and
#'   \code{counts_by_group} is supplied, existing scores are recomputed and
#'   replaced.
#'
#' @return \code{runCCI.blisa}: the input \code{blisa} object with
#'   \code{CCI_scores} populated (a wide data frame -- rows are
#'   \code{"Sender->Receiver"} group pairs, columns are LR pairs).
#'
#' @export
runCCI.blisa <- function(x, counts_by_group = NULL, overwrite = FALSE, ...) {
  if (!is.null(x$CCI_scores)) {
    if (!overwrite) {
      message("CCI_scores already present \u2014 returning object unchanged. ",
              "Set overwrite = TRUE to recompute.")
      return(x)
    }
    if (is.null(counts_by_group))
      stop(
        "counts_by_group is required to recompute CCI_scores. ",
        "Supply counts_by_group explicitly or set overwrite = FALSE to keep existing scores."
      )
    message("Recomputing CCI_scores (overwrite = TRUE)...")
    x$CCI_scores <- runCCI.default(x, counts_by_group)
    message("CCI_scores replaced.")
    return(x)
  }
  if (is.null(counts_by_group))
    stop(
      "counts_by_group is required when CCI_scores have not been computed. ",
      "Re-run blisa() with a 'group' argument, or supply counts_by_group explicitly."
    )
  message("Computing CCI_scores...")
  x$CCI_scores <- runCCI.default(x, counts_by_group)
  message("CCI_scores added.")
  x
}


#' @describeIn runCCI Default method. Performs the raw CCI computation and
#'   returns only the scores data frame. Typically called internally; use
#'   \code{runCCI.blisa} to compute and attach scores to a \code{blisa} object
#'   in one step.
#'
#' @return \code{runCCI.default}: a data frame with \code{"Sender->Receiver"}
#'   row names and one column per significant LR pair containing the
#'   interaction score \code{0.5 * log2(receiver * sender + 1)}.
#'
#' @export
runCCI.default <- function(x, counts_by_group, ...) {
  LRI_sum       <- x$LR_results
  sw            <- x$spatial_weights
  queen_nb_full <- sw$queen_nb_full
  dist_nb_full  <- sw$dist_nb_full
  ct_names      <- names(counts_by_group)

  scores_list <- lapply(seq_len(nrow(LRI_sum)), function(idx) {
    if (LRI_sum$sig_numbers[idx] == 0) return(NULL)

    l_gene <- LRI_sum$ligand.symbol[idx]
    r_gene <- LRI_sum$receptor.symbol[idx]
    mode   <- LRI_sum$ccc_mode[idx]
    sigHH  <- LRI_sum$sig_index[[idx]]

    nb_list  <- if (mode == "nearby") queen_nb_full else dist_nb_full
    sigHH_ng <- unique(c(sigHH, unlist(nb_list[sigHH])))

    # Sum receptor counts in hotspot bins per group (receivers)
    receiver_sums <- sapply(ct_names, function(ct) {
      sum(get_min_expr(r_gene, counts_by_group[[ct]])[sigHH])
    })

    # Sum ligand counts in hotspot + neighbour bins per group (senders)
    sender_sums <- sapply(ct_names, function(ct) {
      sum(get_min_expr(l_gene, counts_by_group[[ct]])[sigHH_ng])
    })

    score_mat  <- 0.5 * log2(outer(receiver_sums, sender_sums, FUN = "*") + 1)
    pair_names <- as.vector(outer(ct_names, ct_names, function(r, s) paste(s, r, sep = "->")))
    scores     <- as.vector(score_mat)
    names(scores) <- pair_names
    scores
  })
  names(scores_list) <- rownames(LRI_sum)

  scores_list <- Filter(Negate(is.null), scores_list)
  if (length(scores_list) == 0) return(data.frame())

  score_df <- as.data.frame(do.call(cbind, scores_list))
  # Split "Sender->Receiver" row names into explicit columns
  sr_pairs <- strsplit(rownames(score_df), "->")
  data.frame(
    Sender   = sapply(sr_pairs, `[`, 1),
    Receiver = sapply(sr_pairs, `[`, 2),
    score_df,
    row.names = NULL,
    check.names = FALSE
  )
}
