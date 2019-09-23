#' Calculate normalized counts.
#'
#' Calculate normalized counts.
#'
#' Normalization uses methods from the package corresponding to the de_pkg
#' paramter:
#' \describe{ \item{edger}{Uses edgeR::cpm with normalized.lib.sizes =
#' TRUE and log = TRUE}
#' \item{deseq2}{Not implemented yet.} }
#'
#' @param x A BbcSE object or a DGEList or a DESeqDataSet.
#' @param de_pkg "edger" or "deseq2". Only used if x is a BbcSE
#' @param ... Not used currently
#' @return A BbcSE object or a matrix or....
#' @seealso \code{\link[edgeR]{cpm}}
#' @export
setGeneric("normalize_counts", function(x, ...)
  standardGeneric("normalize_counts")
)
