#' BbcSE, an extension of SummarizedExperiment for VAI BBC RNA-seq workflow
#'
#' In a \code{BbcSE} object, "counts" must be the first assay and must contain
#' non-negative values.
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SummarizedExperiment assayNames assay
.BbcSE <- setClass("BbcSE", contains="SummarizedExperiment")

setValidity("BbcSE", function(object) {
  msg <- NULL

  if (assayNames(object)[1] != "counts") {
    msg <- c(msg, "'counts' must be first assay")
  }

  if (min(assay(object)) < 0) {
    msg <- c(msg, "negative values in 'counts'")
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

### ---------------------------------------------------------------------

#' BbcRNAData, an S4 class to represent RNA-seq data in various containers
#'
#' @slot summexp A RangedSummarizedExperiment.
#' @slot edger A list containing DGEList (must be first element) and edgeR
#'   results objects.
#' @slot deseq2 A DESeqDataSet object.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importClassesFrom edgeR DGEList DGEGLM DGEExact DGELRT
.BbcRNAData <- setClass("BbcRNAData",
                    slots = list(bbcse = "BbcSE",
                                 edger = "list",
                                 deseq2 = "DESeqDataSet")
)

setValidity("BbcRNAData", function(object) {
  msg <- NULL

  if (length(object@edger) > 0 & object@edger[1] != "DGEList") {
    msg <- c(msg, "'DGEList' must be first element in 'edger' slot.")
  }

  if (length(object@edger) > 1) {
    for (i in 2:length(object@edger)){
      if (!object@edger[i] %in% c("DGEGLM", "DGEExact", "DGELRT")) {
        msg <- c(msg,
                 "After 'DGEList', 'edger' slot must be EdgeR result objects.")
      }
    }
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

#' Constructor for BbcRNAData
#'
#' \code{BbcRNAData} is a constructor for BbcRNAData.
#'
#' @param counts A count matrix with sample names as column names and gene names
#'   as row names.
#' @param ... Arguments for SummarizedExperiment.
#' @export
BbcRNAData <- function(counts, ...) {
  .BbcRNAData(bbcse = BbcSE(counts, ...))
}

### ---------------------------------------------------------------------



