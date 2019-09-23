###-----------------------------------------------------------------------------
#' Getter for edger slot
#' @param x BbcSE object
#' @export
edger <- function(x) {
  if(!is(x, "BbcSE")) stop("Not BbcSE object")
  out <- x@edger
  out
}

#' Setter for edger slot
#' @param x BbcSE object
#' @param value a list of edgeR objects
#' @export
`edger<-` <- function(x, value) {
  if(!is(x, "BbcSE")) stop("x not BbcSE object")

  x@edger <- value

  validObject(x)

  x
}

###-----------------------------------------------------------------------------

#' Make a DGEList object.
#'
#' Make a DGEList object based on assay(object, "counts").
#'
#' @param x A BbcSE object.
#' @param group Name of a column from colData()
#' @param rm_low_genes logical indicating whether low count genes should be
#'   removed. Implements edgeR::filterByExpr.
#' @param calc_norm logical indicating whether normalization factors should be
#'   calculated. Implements edgeR::calcNormFactors. Also turns on calculations
#'   of normalized counts and per-group normalized counts
#' @return A BbcSE object.
#' @seealso \code{\link[edgeR]{filterByExpr}}
#' @importFrom edgeR DGEList filterByExpr calcNormFactors cpm
#' @importFrom SummarizedExperiment assay colData SummarizedExperiment
#' @export
makeDGEList <- function(x, group = NULL, rm_low_genes = TRUE,
                        calc_norm = TRUE) {
  if(!is(x, "BbcSE")) stop("x not BbcSE object")
  if(length(edger(x)) > 0) stop("edger slot not empty")
  if(missing(group)) stop("Please provide group parameter")
  if(!group %in% colnames(colData(x))) stop("Provided group not in colData")

  # make DGEList
  myDGEList <- edgeR::DGEList(counts = assay(x, "counts"),
                             group = colData(x)[[group]])

  if (isTRUE(rm_low_genes)){
    # filter out lowly expressed genes
    rows_keep <- edgeR::filterByExpr(myDGEList,
                                     group = myDGEList$samples$group)

    myDGEList <- myDGEList[rows_keep, ]
  }

  if (isTRUE(calc_norm)){
    # calculate norm factors
    myDGEList <- edgeR::calcNormFactors(myDGEList)

  }

  # store DGEList in the BbcSE object
  edger(x) <- list(DGEList = myDGEList)

  # add log normalized counts if norm factors calculated
  if (isTRUE(calc_norm)) {
    x <- normalize_counts(x)
  }

  validObject(x)

  return(x)
}

###-----------------------------------------------------------------------------
#' @describeIn normalize_counts For de_pkg="edger", a SummarizedExperiment
#'   object is created and stored in an element named "norm_cts" in the edger
#'   slot. Group average log normalized counts are also calculated and stored as
#'   rowData.
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("normalize_counts", "BbcSE", function(x, de_pkg = "edger") {

  if(!de_pkg %in% c("edger", "deseq2")){
    stop("de_pkg not recognized")
  }

  if(identical(de_pkg, "edger")){
    if(length(edger(x)) == 0){
      stop("No DGEList in edger slot. Run makeDGEList first.")
    }

    # calculate normalized counts. Use normalize_counts,DGEList
    norm_log_cpm <- normalize_counts(edger(x)$DGEList)

    # calculate by group cpms
    group_norm_log_cpm <- edgeR::cpmByGroup(edger(x)$DGEList,
                                            normalized.lib.sizes = TRUE,
                                            log=TRUE)

    ## modify col names of group cpms
    colnames(group_norm_log_cpm) <- paste0(colnames(group_norm_log_cpm),
                                           ".norm_log_cpm")

    ## add group cpms to rowData
    new_row_data <- cbind(rowData(x)[rownames(edger(x)$DGEList), , drop = FALSE],
                          group_norm_log_cpm)

    # store normalized counts in a SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(norm_log_cpm = norm_log_cpm),
      rowData = new_row_data,
      colData = colData(x)
    )

  }

  x@edger$norm_cts <- se

  validObject(x)

  return(x)
})

###-----------------------------------------------------------------------------
#' @describeIn normalize_counts log normalized counts returned as a matrix.
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("normalize_counts", "DGEList", function(x) {
  # calculate normalized counts
  norm_log_cpm <- edgeR::cpm(x,
                             normalized.lib.sizes = TRUE,
                             log=TRUE)

  return(norm_log_cpm)
})

