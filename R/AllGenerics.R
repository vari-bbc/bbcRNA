#' Calculate normalized counts.
#'
#' Calculate normalized counts.
#'
#' Normalization uses methods from the package corresponding to the de_pkg
#' paramter:
#' \describe{ \item{edger}{Uses edgeR::cpm with normalized.lib.sizes =
#' TRUE and log = TRUE}
#' \item{deseq2}{Not implemented yet.} }
#'
#' @param x A BbcSE object or a DGEList or a DESeqDataSet.
#' @param de_pkg "edger" or "deseq2". Only used if x is a BbcSE
#' @param ... Not used currently
#' @return A BbcSE object or a matrix or....
#' @seealso \code{\link[edgeR]{cpm}}
#' @export
setGeneric("normalize_counts", function(x, ...)
  standardGeneric("normalize_counts")
)
###-----------------------------------------------------------------------------
#' Identify DE genes.
#'
#' Identify DE genes.
#'
#' Uses methods from the package corresponding to the de_pkg paramter:
#' \describe{ \item{edger}{Uses edgeR::glmQLFit folllowed by either glmQLFTest
#' or glmTreat} \item{deseq2}{Not implemented yet.} }
#'
#' @param x A BbcSE object or a DGEList or a DESeqDataSet.
#' @param de_pkg "edger" or "deseq2". Only used if x is a BbcSE
#' @param design chr value. For example, '~0+group'. Variables in the design
#'   must be present in colData. For de_pkg="edger", passed to glmQLFit
#' @param contrasts list of chr vectors containing variable name, numerator
#'   level, denominator level. For nested contrasts, use the format
#'   'level2-level1' for the numerator and denominator values.
#' @param coefs  list of integers or character vectors indicating which
#'   coefficients of the linear model are to be tested equal to zero. Values
#'   must be columns or column names of design.
#' @param test For de_pkg="edger", either "glmQLFTest" or "glmTreat"
#' @param sample_meta Column meta data as DataFrame or data.frame
#' @param lfc See edgeR::glmTreat. Only used for test="glmTreat".
#' @param ... Not used currently.
#' @return A BbcSE object or a BbcEdgeR object or...
#' @seealso \code{\link[edgeR]{glmQLFit} \link[edgeR]{glmQLFTest}
#'   \link[edgeR]{glmTreat} }
#' @export
setGeneric("findDEGs", function(x, ...)
  standardGeneric("findDEGs")
)

###-----------------------------------------------------------------------------


#' Run GSEA
#'
#' Run GSEA
#'
#'
#' @param x A BbcSE object or a DGEList or a DESeqDataSet.
#' @param de_pkg "edger" or "deseq2". Only used if x is a BbcSE
#' @param gene_set one of "kegg", "reactome", "H", "C1", "C2", "C3", "C4", "C5",
#'   "C6", "C7"
#' @param orgDb OrgDB object for your organism
#' @param organism For compatibility with dowstream tools, the organism is
#'   specified differently depending on the desired gene set. For KEGG, use the
#'   three-letter code according to
#'   http://www.genome.jp/kegg/catalog/org_list.html. For Reactome, possible
#'   values are ‘celegans’, ‘fly’, ‘human’, ‘mouse’, ‘rat’, ‘yeast’ and
#'   ‘zebrafish’. For msigdb, run 'msigdbr::msigdbr_show_species()'.
#' @param contrast_names character value or vector for the contrast(s) of
#'   interest. See name(de_results(edger(x))). If left to default "" value then
#'   all contrasts will be processed.
#' @param rank_by \describe{
#'  \item{"signed-log10pval"}{-log10 of the PValue from the DE analysis, signed
#'  by the LFC.}
#'  \item{"-log10pval"}{-log10 of the PValue from the DE analysis. May be ideal
#'  for gene sets that contain both up-regulated and down-regulated genes for a
#'  particular process or pathway etc.}
#' }
#' @param ... Passed to gseKEGG, GSEA or gsePathway. See each function for
#'   possible arguments.
#' @return A list of gseaResult objects or one gseaResult object
#' @seealso \code{\link[clusterProfiler]{gseKEGG} \link[clusterProfiler]{GSEA}
#' \link[ReactomePA]{gsePathway}}
#' @export
setGeneric("run_gsea", function(x, ...)
  standardGeneric("run_gsea")
)



