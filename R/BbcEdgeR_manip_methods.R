###-----------------------------------------------------------------------------
#' Constructor for BbcEdgeR
#'
#' \code{BbcEdgeR} is a constructor for BbcEdgeR
#'
#' @param dgelist A DGEList
#' @param de_results A list containing DGEGLM as first element and then edgeR
#'   result objects
#' @param norm_cts A SummarizedExperiment
#' @return A BbcEdgeR object.
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom edgeR DGEList
#' @importFrom testthat expect_warning
BbcEdgeR <- function(dgelist = NULL,
                     de_results = list(),
                     norm_cts = SummarizedExperiment())
{
  if (is.null(dgelist)){
    testthat::expect_warning(dgelist <- edgeR::DGEList(),
                             "no non-missing arguments to min; returning Inf")
  }

  bbcedger_obj <- .BbcEdgeR( dgelist = dgelist,
                             de_results = de_results,
                             norm_cts = norm_cts)

  validObject(bbcedger_obj)

  bbcedger_obj
}
###-----------------------------------------------------------------------------

#' Getter for dgelist slot from BbcEdgeR object
#' @param x BbcEdgeR object
#' @export
dgelist <- function(x) {
  if(!is(x, "BbcEdgeR")) stop("x is not a BbcEdgeR object")

  out <- x@dgelist

  out
}

###-----------------------------------------------------------------------------
#' Setter for dgelist slot in a BbcEdgeR object
#'
#' @param x A BbcEdgeR object.
#' @param value A DGEList object
#' @export
`dgelist<-` <- function(x, value) {
  if(!is(x, "BbcEdgeR")) stop("x is not a BbcEdgeR object")

  # set the dgelist slot.
  x@dgelist <- value

  validObject(x)
  x
}

###-----------------------------------------------------------------------------
#' Getter for de_results slot from BbcEdgeR object
#' @param x BbcEdgeR object
#' @export
de_results <- function(x) {
  if(!is(x, "BbcEdgeR")) stop("x is not a BbcEdgeR object")

  out <- x@de_results

  out
}

###-----------------------------------------------------------------------------
#' Setter for de_results slot in a BbcEdgeR object
#'
#' @param x A BbcEdgeR object.
#' @param value A list containing DGEGLM in first element and then edgeR result
#'   objects
#' @export
`de_results<-` <- function(x, value) {
  if(!is(x, "BbcEdgeR")) stop("x is not a BbcEdgeR object")

  # set the de_results slot.
  x@de_results <- value

  validObject(x)
  x
}

###-----------------------------------------------------------------------------
#' Getter for norm_cts slot from BbcEdgeR object
#' @param x BbcEdgeR object
#' @export
norm_cts <- function(x) {
  if(!is(x, "BbcEdgeR")) stop("x is not a BbcEdgeR object")

  out <- x@norm_cts

  out
}

###-----------------------------------------------------------------------------
#' Setter for norm_cts slot in a BbcEdgeR object
#'
#' @param x A BbcEdgeR object.
#' @param value A list containing DGEGLM in first element and then edgeR result
#'   objects
#' @export
`norm_cts<-` <- function(x, value) {
  if(!is(x, "BbcEdgeR")) stop("x is not a BbcEdgeR object")

  # set the norm_cts slot.
  x@norm_cts <- value

  validObject(x)
  x
}
