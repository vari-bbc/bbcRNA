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

#' Make a DGEList object.
#'
#' Make a DGEList object based on assay(object, "counts").
#'
#' @param x A BbcSE object.
#' @param group Name of a column from colData()
#' @param rm_low_genes logical indicating whether low count genes should be
#'   removed. Implements edgeR::filterByExpr.
#' @param calc_norm logical indicating whether normalization factors should be
#'   calculated. Implements edgeR::calcNormFactors.
#' @return A BbcSE object.
#' @seealso \code{\link[edgeR]{filterByExpr}}
#' @importFrom edgeR DGEList filterByExpr calcNormFactors
#' @importFrom SummarizedExperiment assay colData
#' @export
makeDGEList <- function(x, group = NULL, rm_low_genes = TRUE,
                        calc_norm = TRUE) {
  if(!is(x, "BbcSE")) stop("x not BbcSE object")
  if(missing(group)) stop("Please provide group parameter")

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

  # store final DGEList in edger()[[1]] of the BbcSE object
  edger(x) <- list(myDGEList)

  validObject(x)

  return(x)
}

