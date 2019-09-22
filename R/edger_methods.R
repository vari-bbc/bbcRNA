#' Getter for edger slot
#' @name edger
#' @param x BbcSE object
#' @export
edger <- function(x) {
  if(!is(x, "BbcSE")) stop("Not BbcSE object")
  out <- x@edger
  out
}

#' Make a DGEList object.
#'
#' Make a DGEList object based on assay(object, "counts").
#'
#' @name makeDGEList
#' @param x A BbcSE object.
#' @param group Name of a column from colData()
#' @return A BbcSE object.
#' @importFrom edgeR DGEList
#' @importFrom SummarizedExperiment assay colData
#' @export
makeDGEList <- function(x, group = NULL) {
  metadata(x)$edger[[1]] <- edgeR::DGEList(
    counts = assay(x, "counts"),
    group = colData(x)[[group]])

  return(x)
}

