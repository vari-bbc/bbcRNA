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
#' @param aln_rates Matrix containing alignment metrics. Optional.
#' @param ... Not used currently.
#' @return A BbcSE object (extension of RangedSummarizedExperiment).
#' @export
#' @importFrom methods as is
#' @importFrom GenomicRanges GRanges
#' @importFrom stringr str_remove
#' @importFrom SummarizedExperiment SummarizedExperiment
BbcSE <- function(counts = matrix(0, 0, 0),
                  granges = rep(GRanges(seqnames="foobar", ranges=0:0),
                                nrow(counts)),
                  aln_rates = matrix(0, 0, 0),
                  ...)
{
  stopifnot(is(counts, "matrix"))

  if(!is.null(names(granges))){
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
  } else {
    # default granges is fake GRanges data. Needed to get DEFormats::DGEList to
    # work
    se <- SummarizedExperiment(list(counts=counts),
                               rowRanges = granges)
  }

  if (length(aln_rates) > 0){
    # check that sample names match
    stopifnot(identical(nrow(aln_rates), ncol(counts)))
    stopifnot(all(rownames(aln_rates) %in% colnames(counts)))

    # arrange aln_rates rows into same order as counts
    aln_rates <- aln_rates[colnames(counts), ]
  }

  bbcse_obj <- .BbcSE(se,
                      aln_rates = aln_rates,
                      edger = BbcEdgeR(),
                      deseq2 = list())

  validObject(bbcse_obj)

  bbcse_obj
}

###-----------------------------------------------------------------------------
#' show method for BbcSE
#' @param object BbcSE object
#' @export
#' @importMethodsFrom SummarizedExperiment show
#' @importFrom methods slotNames
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

  bbcedger_obj <- edger(object)
  bbcedger_obj_slots <- methods::slotNames(bbcedger_obj)
  cat(paste0("edger(", length(bbcedger_obj_slots), "):"), bbcedger_obj_slots,"\n")

  deseq2_classes <- sapply(deseq2(object), class)
  cat(paste0("deseq2(", length(deseq2_classes), "):"), deseq2_classes,"\n")
})

###-----------------------------------------------------------------------------
#' Subsetting method for BbcSE
#' @param x BbcSE object
#' @param i row index
#' @param j column index
#' @param drop see help("[")
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
  BiocGenerics:::replaceSlots(out, aln_rates = aln_rates,
                              edger = list(),
                              deseq2 = list(),
                              check = TRUE)
})

###-----------------------------------------------------------------------------
#' Assigning subsets method for BbcSE
#' @param x BbcSE object
#' @param i row index
#' @param j column index
#' @param ... see ?"[" for S4Vectors
#' @param value BbcSE object for replacement
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
                                               edger = list(),
                                               deseq2 = list(),
                                               check=TRUE)
                 })

###-----------------------------------------------------------------------------
#' Row combining method for BbcSE
#' @param ... BbcSE objects
#' @param deparse.level see ?`rbind,SummarizedExperiment-method`
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
  BiocGenerics:::replaceSlots(out, edger = list(), deseq2 = list(), check=FALSE)
})

###-----------------------------------------------------------------------------
#' Column combining method for BbcSE
#' @param ... BbcSE objects
#' @param deparse.level See ?base::cbind for a description of this argument.
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
                              edger = list(),
                              deseq2 = list(),
                              check=FALSE)
})
