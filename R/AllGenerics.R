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
###-----------------------------------------------------------------------------


#' Identify DE genes.
#'
#' Identify DE genes.
#'
#' Uses methods from the package corresponding to the de_pkg paramter:
#' \describe{ \item{edger}{Uses edgeR::glmQLFit folllowed by either glmQLFTest
#' or glmTreat} \item{deseq2}{Not implemented yet.} }
#'
#' @param x A BbcSE object or a DGEList or a DESeqDataSet.
#' @param de_pkg "edger" or "deseq2". Only used if x is a BbcSE
#' @param design chr value. For example, '~0+group'. Variables in the design
#'   must be present in colData. For de_pkg="edger", passed to glmQLFit
#' @param contrasts list of chr vectors containing variable name, numerator
#'     level, denominator level
#' @param test For de_pkg="edger", either "glmQLFTest" or "glmTreat"
#' @param sample_meta Column meta data as DataFrame or data.frame
#' @param lfc See edgeR::glmTreat. Only used for test="glmTreat".
#' @return A BbcSE object or a list of DGElist and edgeR result objects or...
#' @seealso \code{\link[edgeR]{glmQLFit} \link[edgeR]{glmQLFTest}
#'   \link[edgeR]{glmTreat} }
#' @export
setGeneric("findDEGs", function(x, ...)
  standardGeneric("findDEGs")
)

