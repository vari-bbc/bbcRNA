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
.BbcSE <- setClass("BbcSE",
                   slots = representation(
                     aln_rates = "matrix"),
                   contains="RangedSummarizedExperiment")

setValidity("BbcSE", function(object) {
  NR <- nrow(object)
  NC <- ncol(object)
  msg <- NULL

  if (assayNames(object)[1] != "counts") {
    msg <- c(msg, "'counts' must be first assay")
  }

  if (min(assay(object)) < 0) {
    msg <- c(msg, "negative values in 'counts'")
  }

  if (!identical(names(metadata(object))[1:2], c("edger", "deseq2"))) {
    msg <- c(msg, "metadata(object)[1:2] names must be 'edger' and 'deseq2'")
  }

  if (!is(metadata(object)[[1]], "list")) {
    msg <- c(msg, "First element in metadata must be a list of edgeR objects")
  }

  if (length(metadata(object)[[1]]) > 0 &&
      !is(metadata(object)[[1]][[1]], "DGEList")) {
    msg <- c(msg, "metadata(object)[[1]][[1]] must be a DGEList object")
  }

  if (length(metadata(object)[[1]]) > 1) {
    for (i in 2:length(metadata(object)[[1]])){
      curr_meta <- metadata(object)[[1]][[i]]
      if (!is(curr_meta, "DGEGLM") &
          !is(curr_meta, "DGEExact") &
          !is(curr_meta, "DGELRT")) {
        msg <- c(msg,
                 "After 'DGEList', metadata(object)[[1]] elements must be
                 edgeR result objects.")
      }
    }
  }

  if (length(metadata(object)) > 1){
    if (!is(metadata(object)[[2]], "list")) {
      msg <- c(msg, "Second element in metadata must be a list of DESeq2 objects")
    }

    if (length(metadata(object)[[2]]) > 0 &&
        !is(metadata(object)[[2]][[1]], "DESeqDataSet")) {
      msg <- c(msg, "metadata(object)[[2]][[1]] must be a DESeqDataSet object")
    }

    if (length(metadata(object)[[2]]) > 1 &&
        !is(metadata(object)[[2]][[2]], "DESeqResults")) {
      msg <- c(msg, "metadata(object)[[2]][[2]] must be a DESeqResults object")
    }
  }


  if (!is.matrix(aln_rates(object, withDimnames=FALSE))) {
    msg <- c(msg, "aln_rates must be a matrix")
  }

  valid_aln_rates_colnames <- c("input_reads",
                                "uniq_aln_reads",
                                "mult_aln_reads",
                                "map_rate",
                                "uniq_map_rate")
  if(all(colnames(aln_rates(object, withDimnames=FALSE)) %in%
         valid_aln_rates_colnames)) {

  }

  if (length(aln_rates(object, withDimnames=FALSE) > 0)) {
    if (nrow(aln_rates(object, withDimnames=FALSE)) != NC) {
      msg <- c(
        msg, "'nrow(aln_rates)' should be equal to the number of columns"
      )
    }

  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

###-----------------------------------------------------------------------------
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
#' @importFrom GenomicRanges GRanges
#' @importFrom stringr str_remove
#' @importFrom SummarizedExperiment SummarizedExperiment
BbcSE <- function(counts = matrix(0, 1, 1),
                  granges,
                  aln_rates = matrix(0, 0, 0),
                  ...)
{
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
                               rowRanges = granges[common_genes],
                               metadata = list(edger = list(), deseq2 = list()))
  }

  if(missing(granges)){
    # add fake GRanges data. Needed to get DEFormats::DGEList to work
    se <- SummarizedExperiment(list(counts=counts),
                               rep(GRanges(seqnames="foobar", ranges=0:0),
                                   nrow(counts)),
                               metadata = list(edger = list(), deseq2 = list()))

    # se <- as(se, "RangedSummarizedExperiment")
  }

  if (length(aln_rates) > 0){
    # check that sample names match
    stopifnot(identical(nrow(aln_rates), ncol(counts)))
    stopifnot(all(rownames(aln_rates) %in% colnames(counts)))

    # arrange aln_rates rows into same order as counts
    aln_rates <- aln_rates[colnames(counts), ]
  }

  .BbcSE(se, aln_rates = aln_rates)
}

###-----------------------------------------------------------------------------
#' show method for BbcSE
#'
#' @export
#' @importMethodsFrom SummarizedExperiment show
setMethod("show", "BbcSE", function(object) {
  out <- callNextMethod()

  aln_rates_ncol <- ncol(aln_rates(object, withDimnames=FALSE))
  aln_rates_colnames <- ifelse(aln_rates_ncol > 0,
                               paste(colnames(
                                 aln_rates(object, withDimnames=TRUE)
                               ), collapse = " "),
                               "")

  cat(
    "aln_rates(", aln_rates_ncol, "): ", aln_rates_colnames, "\n",
    sep=""
  )
})

###-----------------------------------------------------------------------------
#' Subsetting method for BbcSE
#'
#' @export
setMethod("[", "BbcSE", function(x, i, j, drop=TRUE) {
  aln_rates <- aln_rates(x, withDimnames=FALSE)

  if (!missing(i)) {
    if (is.character(i)) {
      fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
      i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        i, rownames(x), fmt
      )
    }
    i <- as.vector(i)
  }

  if (!missing(j)) {
    if (is.character(j)) {
      fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
      j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
        j, colnames(x), fmt
      )
    }
    j <- as.vector(j)
    aln_rates <- aln_rates[j, , drop=FALSE]
  }

  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, aln_rates = aln_rates, check = TRUE)
})

###-----------------------------------------------------------------------------
#' Assigning subsets method for BbcSE
#' @export
setReplaceMethod("[", c("BbcSE", "ANY", "ANY", "BbcSE"),
                 function(x, i, j, ..., value) {

                   aln_rates <- aln_rates(x, withDimnames=FALSE)

                   if (!missing(i)) {
                     if (is.character(i)) {
                       fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
                       i <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                         i, rownames(x), fmt
                       )
                     }
                     i <- as.vector(i)

                   }

                   if (!missing(j)) {
                     if (is.character(j)) {
                       fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
                       j <- SummarizedExperiment:::.SummarizedExperiment.charbound(
                         j, colnames(x), fmt
                       )
                     }
                     j <- as.vector(j)

                     aln_rates[j,] <- aln_rates(value, withDimnames=FALSE)
                   }

                   out <- callNextMethod()
                   BiocGenerics:::replaceSlots(out,
                                               aln_rates=aln_rates,
                                               check=TRUE)
                 })

###-----------------------------------------------------------------------------
#' Row combining method for BbcSE
#' @importMethodsFrom SummarizedExperiment rbind
#' @export
setMethod("rbind", "BbcSE", function(..., deparse.level=1) {
  args <- list(...)

  # Checks for identical column state.
  ref <- args[[1]]
  ref.aln_rates <- aln_rates(ref, withDimnames=FALSE)
  for (x in args[-1]) {
    if (!identical(ref.aln_rates, aln_rates(x, withDimnames=FALSE)))
    {
      stop("per-column values are not compatible")
    }
  }

  old.validity <- S4Vectors:::disableValidity()
  S4Vectors:::disableValidity(TRUE)
  on.exit(S4Vectors:::disableValidity(old.validity))

  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, check=FALSE)
})

###-----------------------------------------------------------------------------
#' Column combining method for BbcSE
#' @importMethodsFrom SummarizedExperiment rbind
#' @export
setMethod("cbind", "BbcSE", function(..., deparse.level=1) {
  args <- list(...)

  all.aln_rates <- lapply(args, aln_rates, withDimnames=FALSE)

  all.aln_rates <- do.call(rbind, all.aln_rates)

  # Checks for identical column state.
  # ref <- args[[1]]
  # ref.rv <- rowVec(ref, withDimnames=FALSE)
  # ref.rrm <- rowToRowMat(ref, withDimnames=FALSE)
  # ref.crm <- colToRowMat(ref, withDimnames=FALSE)
  # for (x in args[-1]) {
  #   if (!identical(ref.rv, rowVec(x, withDimnames=FALSE))
  #       || !identical(ref.rrm, rowToRowMat(x, withDimnames=FALSE))
  #       || !identical(ref.crm, colToRowMat(x, withDimnames=FALSE)))
  #   {
  #     stop("per-row values are not compatible")
  #   }
  # }

  old.validity <- S4Vectors:::disableValidity()
  S4Vectors:::disableValidity(TRUE)
  on.exit(S4Vectors:::disableValidity(old.validity))

  out <- callNextMethod()
  BiocGenerics:::replaceSlots(out, aln_rates=all.aln_rates,
                              check=FALSE)
})
