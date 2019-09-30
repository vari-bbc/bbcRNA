#' Prep clusterProfiler geneList
#'
#' Prep clusterProfiler geneList
#' Adapts code from
#' https://github.com/YuLab-SMU/DOSE/wiki/how-to-prepare-your-own-geneList.
#'
#'
#' @param gene_and_metric_df two column df with Entrez ID as col1 and
#'   metric(log2FC or F etc) as col2
#' @export
prep_clusterprofiler_genelist <- function(gene_and_metric_df){
  ## feature 1: numeric vector
  geneList <- gene_and_metric_df[[2]]
  ## feature 2: named vector
  names(geneList) <- as.character(gene_and_metric_df[[1]])
  ## feature 3: decreasing order
  geneList <- sort(geneList, decreasing = TRUE)

  return(geneList)
}

###-----------------------------------------------------------------------------
#' @describeIn run_gsea x is a dataframe with gene name as col1 and
#'   metric(log2FC or F etc) as col2. Function will sort by col2 in decreasing
#'   order. See '?prep_clusterprofiler_genelist'.
#' @importFrom clusterProfiler GSEA gseKEGG
#' @importFrom ReactomePA gsePathway
#' @importFrom DOSE setReadable
#' @importFrom msigdbr msigdbr
#' @export
setMethod("run_gsea", "data.frame", function(x, gene_set="reactome", organism,
                                         orgDb, ...) {

  ellipses_forbidden_args <- c("geneList", "keyType", "TERM2GENE", "TERM2NAME")
  args <- list(...)
  args_forbidden <- args %in% ellipses_forbidden_args
  if(any(args_forbidden)) {
    stop(paste0("args forbidden: ", paste0(unlist(args[args_forbidden]),
                                           collapse = " ")
    ))
  }

  # prep geneList
  gene_list <- prep_clusterprofiler_genelist(x)

  if (identical("reactome", gene_set)){

    gsea_res <- ReactomePA::gsePathway(geneList = gene_list,
                                       organism = organism,
                                       ...)

  } else if (identical("kegg", gene_set)){

    gsea_res <- clusterProfiler::gseKEGG(geneList = gene_list,
                                         organism = organism,
                                         keyType = "kegg",
                                         ...)

  } else if (gene_set %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")){
    msigdbr_cat <- gene_set
    msigdbr_gene_set <- msigdbr::msigdbr(species = organism,
                                         category = gene_set)

    # convert to clusterprofiler friendly format
    msigdbr_clusterprofiler <- msigdbr_gene_set[, c("gs_name", "entrez_gene")] %>%
      as.data.frame(stringsAsFactors = FALSE)

    gsea_res <- clusterProfiler::GSEA(geneList = gene_list,
                                      TERM2GENE = msigdbr_clusterprofiler,
                                      ...)

  } else{
    stop("Invalid value for gene_set parameter. See ?run_gsea.")
  }

  gsea_res <- DOSE::setReadable(gsea_res,
                                OrgDb = orgDb,
                                keyType = "ENTREZID")

  gsea_res
})

###-----------------------------------------------------------------------------
#' @describeIn run_gsea Extract out the gene and signed -log10 PValue. Run
#'   run_gsea for each result object (contrast).
#' @export
setMethod("run_gsea", "BbcSE", function(x, de_pkg="edger",
                                        contrast_names = "", ...) {
  if(!"entrez_ids" %in% colnames(rowData(x))) {
    stop("Please run ens2entrez first to get 'entrez_ids' column in rowData")
  }

  if(identical("edger", de_pkg)){
    # get all the contrasts
    edger_results <- de_results(edger(x))[-1] # first element is a DGEGLM object

    if(identical(contrast_names, "")){
      edger_results_names <- names(edger_results)
    } else {
      if(!all(contrast_names %in% names(edger_results))){
        stop("Not all contrast names found.")
      }
      # get requested contrasts
      edger_results <- de_results(edger(x))[contrast_names]
      edger_results_names <- contrast_names
    }

    gsea_results <- lapply(edger_results_names, function(edger_res_name){
      dge_table <- edger_results[[edger_res_name]]$table
      dge_table$entrez_ids <- rowData(x)$entrez_ids[rownames(dge_table)]

      # remove genes with logFC == 0 or PValue == 1. Pvalue = 1 work out to
      # minuslog10pval of 0, so not very useful anyways. logFC == 0 cannot be
      # considered up-regulated or down-regulated
      rows_b4_filt <- nrow(dge_table)
      filt_dge_table <- dge_table[dge_table$logFC != 0 &
                                    dge_table$PValue != 1, ]
      message(paste0(edger_res_name, ": removed ",
                     rows_b4_filt - nrow(filt_dge_table),
                     " genes out of ", rows_b4_filt,
                     " due to logFC = 0 or PValue = 1"))

      # remove genes with no Entrez match (is.na)
      rows_b4_filt <- nrow(filt_dge_table)
      filt_dge_table <- filt_dge_table[!is.na(filt_dge_table$entrez_ids), ]
      diff <- rows_b4_filt - nrow(filt_dge_table)
      message(paste0(edger_res_name, ": removed ",
                     rows_b4_filt - nrow(filt_dge_table),
                     " genes out of ", rows_b4_filt,
                     " due to no Entrez gene match"))

      message(paste0(edger_res_name, ": total genes remaining is ",
                     nrow(filt_dge_table)))

      # calculate -log10(Pvalue) signed by LFC
      filt_dge_table$signed_minuslog10Pval <-
        sign(filt_dge_table$logFC) * -log10(filt_dge_table$PValue)

      # get entrezgene and signed_minuslog10Pval
      genes_and_score <- filt_dge_table[,c("entrez_ids",
                                           "signed_minuslog10Pval")]

      # run the data.frame method of 'run_gsea'
      run_gsea(x = genes_and_score, ...)
    })

    names(gsea_results) <- edger_results_names

  }

  return(gsea_results)
})

###-----------------------------------------------------------------------------
