#' Getter for aln_rates slot
#' @export
setMethod("aln_rates", "BbcSE", function(x, withDimnames=TRUE) {
  out <- x@aln_rates
  if (withDimnames) {
    rownames(out) <- colnames(x)
  }
  out
})

#' Setter for aln_rates slot
#'
#' @param x A BbcSE object.
#' @param value A matrix of mapping metrics. Rownames must correspond to samples
#' in counts matrix and column names indicating metric type must be set.
#' @export
setReplaceMethod("aln_rates", "BbcSE", function(x, value) {

    # check that sample names match
  stopifnot(identical(nrow(value), ncol(x)))
  stopifnot(all(rownames(value) %in% colnames(x)))

  # set the aln_rates slot. Order the aln_rates according to colnames(x)
  x@aln_rates <- value[colnames(x), ]
  validObject(x)
  x
})
