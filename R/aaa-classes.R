#' BbcSE, an extension of SummarizedExperiment for VAI BBC RNA-seq workflow
#'
#' In a \code{BbcSE} object, "counts" must be the first assay and must contain
#' non-negative values.
#'
#' @importFrom S4Vectors metadata
#' @importFrom  SummarizedExperiment assayNames assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importClassesFrom edgeR DGEList DGEGLM DGEExact DGELRT
.BbcSE <- setClass("BbcSE", contains="SummarizedExperiment")

setValidity("BbcSE", function(object) {
  msg <- NULL

  if (assayNames(object)[1] != "counts") {
    msg <- c(msg, "'counts' must be first assay")
  }

  if (min(assay(object)) < 0) {
    msg <- c(msg, "negative values in 'counts'")
  }

  if (length(metadata(object)) > 0){
    if (class(metadata(object)[[1]]) !=  "list") {
      msg <- c(msg, "First element in metadata must be a list of edgeR objects")
    }

    if (length(metadata(object)[[1]]) > 0 &
        class(metadata(object)[[1]][[1]]) !=  "DGEList") {
      msg <- c(msg, "metadata(object)[[1]][[1]] must be a DGEList object")
    }

    if (length(metadata(object)[[1]]) > 1) {
      for (i in 2:length(metadata(object)[[1]])){
        if (!class(metadata(object)[[1]][i]) %in%
            c("DGEGLM", "DGEExact", "DGELRT")) {
          msg <- c(msg,
                   "After 'DGEList', metadata(object)[[1]] elements must be
                 edgeR result objects.")
        }
      }
    }

    if (length(metadata(object)) > 1 &
        class(metadata(object)[[2]]) !=  "DESeqDataSet") {
      msg <- c(msg, "Second element in metadata must be a DESeqDataSet object")
    }
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

#' Constructor for BbcSE
#'
#' \code{BbcSE} is a constructor for BbcSE.
#'
#' @param counts A count matrix with sample names as column names and gene names
#'   as row names.
#' @param ... Arguments for SummarizedExperiment.
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
BbcSE <- function(counts, ...) {
  se <- SummarizedExperiment(list(counts=counts), ...)
  .BbcSE(se)
}

