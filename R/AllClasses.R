#' BbcSE, an extension of SummarizedExperiment for VAI BBC RNA-seq workflow
#'
#' In a \code{BbcSE} object, "counts" must be the first assay and must contain
#' non-negative values.
#' @export
#' @importFrom S4Vectors metadata
#' @importFrom  SummarizedExperiment assayNames assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importClassesFrom edgeR DGEList DGEGLM DGEExact DGELRT
.BbcSE <- setClass("BbcSE", contains="RangedSummarizedExperiment")

setValidity("BbcSE", function(object) {
  msg <- NULL

  if (assayNames(object)[1] != "counts") {
    msg <- c(msg, "'counts' must be first assay")
  }

  if (min(assay(object)) < 0) {
    msg <- c(msg, "negative values in 'counts'")
  }

  if (length(metadata(object)) > 0){
    if (is(metadata(object)[[1]], "list")) {
      msg <- c(msg, "First element in metadata must be a list of edgeR objects")
    }

    if (length(metadata(object)[[1]]) > 0 &
        is(metadata(object)[[1]][[1]], "DGEList")) {
      msg <- c(msg, "metadata(object)[[1]][[1]] must be a DGEList object")
    }

    if (length(metadata(object)[[1]]) > 1) {
      for (i in 2:length(metadata(object)[[1]])){
        curr_meta <- metadata(object)[[1]][i]
        if (!is(curr_meta, "DGEGLM") &
            !is(curr_meta, "DGEExact") &
            !is(curr_meta, "DGELRT")) {
          msg <- c(msg,
                   "After 'DGEList', metadata(object)[[1]] elements must be
                 edgeR result objects.")
        }
      }
    }

    if (length(metadata(object)) > 1 &
        !is(metadata(object)[[2]], "DESeqDataSet")) {
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
#' If GRanges object is provided, the names must match the rownames of the
#'     counts matrix. The output object will show only the intersected genes
#'     between the the counts and the GRanges.
#'
#' @param counts A count matrix with sample names as column names and gene names
#'   as row names.
#' @param granges GRanges or GRangesList object. Optional.
#' @param ... Arguments for SummarizedExperiment.
#' @return A BbcSE object (extension of RangedSummarizedExperiment).
#' @export
#' @importFrom methods as is
#' @import GenomicFeatures
#' @importFrom GenomicRanges GRanges
#' @importFrom stringr str_remove
#' @importFrom SummarizedExperiment SummarizedExperiment
BbcSE <- function(counts, granges, ...) {
  # se <- as(SummarizedExperiment(list(counts=counts), ...),
  #          "RangedSummarizedExperiment")
  # txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file,
  #                                             format = "gtf")
  # granges <- GenomicFeatures::transcriptsBy(txdb, by="gene")
  stopifnot(is(counts, "matrix"))

  if(!missing(granges)){
    stopifnot(is(granges, "GRanges"))
    names(granges) <- stringr::str_remove(names(granges), "\\.\\d+$")
    common_genes <- intersect(rownames(counts), names(granges))

    num_genes <- sapply(list(common_genes, rownames(counts), names(granges)),
                        length)

    if (num_genes[1] != num_genes[2] | num_genes[1] != num_genes[3]){
      warning("Genes in counts and GRanges object not matched. Taking intersect.")
    }

    se <- SummarizedExperiment(list(counts = counts[common_genes, ]),
                               rowRanges = granges[common_genes])
  }

  if(missing(granges)){
    # add fake GRanges data. Needed to get DEFormats::DGEList to work
    se <- SummarizedExperiment(list(counts=counts),
                               rep(GRanges(seqnames="foobar", ranges=0:0),
                                   nrow(counts)))

    # se <- as(se, "RangedSummarizedExperiment")
  }

  return(.BbcSE(se))
}

