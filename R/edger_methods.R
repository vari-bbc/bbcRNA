#' Convert Ensembl IDs to gene symbols.
#'
#' Convert Ensembl IDs to gene symbols.
#'
#' Non 1:1 matches are determined from all the Ensembl IDs in the OrgDb, not
#' just those in the BbcSE object.
#'
#' @name makeDGEList
#' @param x A BbcSE object.
#' @param group Name of a column from colData()
#' @importFrom edgeR DGEList
#' @importFrom SummarizedExperiment assay colData
#' @return A BbcSE object.
setMethod("makeDGEList", c(x = "BbcSE"),
          function(x, group = NULL) {
            metadata(x)$edger[[1]] <- edgeR::DGEList(
              counts = assay(x, "counts"),
              group = colData(x)[[group]])

            return(x)
          }
)
