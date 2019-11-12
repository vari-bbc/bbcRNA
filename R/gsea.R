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
                                        rank_by="signed-log10pval",
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

      # calculate rank_metric
      if(identical(rank_by, "signed-log10pval")){
        filt_dge_table$rank_metric <-
          sign(filt_dge_table$logFC) * -log10(filt_dge_table$PValue)
      } else if(identical(rank_by, "-log10pval")){
        filt_dge_table$rank_metric <- -log10(filt_dge_table$PValue)
      } else{
        stop("Specify valid ranking metric.")
      }

      # get entrezgene and rank_metric
      genes_and_score <- filt_dge_table[,c("entrez_ids",
                                           "rank_metric")]

      # run the data.frame method of 'run_gsea'
      run_gsea(x = genes_and_score, ...)
    })

    names(gsea_results) <- edger_results_names

  }

  return(gsea_results)
})

###-----------------------------------------------------------------------------
#' Compare genes in clusterProfiler gseaResult to genes in a bbcRNA object
#'
#' Output dataframe for genes that did not have an Entrez match (based on rules
#' in ens2entrez) and genes that are present in the bbcRNA object but not in the
#' bbcRNA object. Extracts LFC and PValue from DE results. If present in the
#' bbcRNA object but not in the DE analysis, then it is inferred to be a lowly
#' expressed gene. Includes all rowData.
#'
#' @param gseaResults A list of ClusterProfiler::gseaResult objects outputted by
#'   'run_gsea'.
#' @param bbcRNA  bbcRNA object
#' @param de_pkg "edger" or "deseq2"
#' @return A dataframe of missing genes.
#' @export
#' @importFrom dplyr filter mutate left_join
#' @importFrom tibble rownames_to_column
find_missing_in_gseaResult <- function(gseaResults, bbcRNA, de_pkg="edger") {

  # row data from bbcRNA object
  bbcRNA_rowdata <- rowData(bbcRNA)

  # bbcRNA object must have 'entrez_ids' column
  if(!"entrez_ids" %in% colnames(bbcRNA_rowdata)){
    stop("Must have run ens2entrez with bbcRNA object. Use bbcRNA object used for 'run_gsea'.")
  }

  # compare gsea results to edgeR results
  if(identical(de_pkg, "edger")){
    gseaResults_names <- names(gseaResults)
    names(gseaResults_names) <- gseaResults_names

    # list of contrasts and list of gsea results must have same element names
    if(!all(gseaResults_names %in% names(bbcRNA@edger@de_results))){
      stop("Not all names of gseaResults matched by contrasts in bbcRNA object.")
    }

    # for each gseaResult name, identify the missing genes and get LFC and PValue from corresponding edgeR result
    missing_genes_list <- lapply(gseaResults_names, function(gseaRes_name){

      # subset for the current gseaResult object from the list and get the gene names (which are in Entrez format)
      curr_gseaRes <- gseaResults[[gseaRes_name]]
      curr_gseaRes_entrez <- names(curr_gseaRes@geneList)

      # check that all the gseaResult genes are in the rowdata of the bbcRNA object
      if(!all(curr_gseaRes_entrez %in% bbcRNA_rowdata$entrez_ids)){
        stop("Not all genes in gseaResult object found in bbcRNA object. Be sure to use same bbcRNA object used for 'run_gsea'.")
      }

      # get genes present in bbcRNA object but not gseaResult
      missing_genes <- bbcRNA_rowdata %>%
        as.data.frame(stringsAsFactors=FALSE) %>%
        tibble::rownames_to_column("ens_gene") %>%
        dplyr::filter(!.data$entrez_ids %in% curr_gseaRes_entrez)

      # get the edgeR results corresponding to the current gseaResult
      edger_de_table <- bbcRNA@edger@de_results[[gseaRes_name]]$table %>%
        tibble::rownames_to_column("ens_gene")
      edger_genes <- edger_de_table$ens_gene

      # sanity check that all genes in the edgeR results are in the bbcRNA object. Should never be triggered.
      if(!all(edger_genes %in% rownames(bbcRNA_rowdata))){
        stop("Not all genes in edgeR object found in bbcRNA object.")
      }

      missing_genes_detail <- missing_genes %>%
        # assume genes absent from edger object but present in bbcRNA to be lowly expressed
        dplyr::mutate(low_expr=ifelse(!missing_genes$ens_gene %in% edger_genes,
                                      TRUE, FALSE)) %>%
        # merge missing genes with edgeR LFC and PValue
        dplyr::left_join(edger_de_table[,c("ens_gene","logFC","PValue"), drop=FALSE], by="ens_gene")

      missing_genes_detail
    })
  }

  return(missing_genes_list)
}

###-----------------------------------------------------------------------------
#' @describeIn shorten_desc Shorten description in enrichResult.
#' @export
#' @importFrom stringr str_trunc
setMethod("shorten_desc", "enrichResult", function(x, max_len=50) {
  x@result$Description <- stringr::str_trunc(x@result$Description, max_len)
  return(x)
})

###-----------------------------------------------------------------------------

#' @describeIn shorten_desc Shorten description in enrichResult.
#' @export
#' @importFrom stringr str_trunc
setMethod("shorten_desc", "gseaResult", function(x, max_len=50) {
  x@result$Description <- stringr::str_trunc(x@result$Description, max_len)
  return(x)
})

###-----------------------------------------------------------------------------

#' @describeIn run_hypergeometric x is a dataframe with gene name (entrez) as
#'   col1, log2FC as col2 and adjusted PValue as col3. Function will filter by
#'   LFC and PValue for DE genes, and optionally will split the DE genes into up
#'   and down-regulated based on LFC. The gene universe are all the genes in the
#'   dataframe. Returns a list of enrichResults for up and down-regulaed genes
#'   or just all the DE genes without regard for direction.
#'
#' @importFrom clusterProfiler enrichKEGG enrichGO enricher
#' @importFrom ReactomePA enrichPathway
#' @importFrom DOSE setReadable
#' @importFrom msigdbr msigdbr
#' @export
setMethod("run_hypergeometric", "data.frame",
          function(x, gene_set="reactome", organism, orgDb,
                   split_by_lfc_dir=TRUE, padj_cutoff=0.05, lfc_cutoff=0, ...) {

  ellipses_forbidden_args <- c("geneList", "keyType", "TERM2GENE", "TERM2NAME")
  args <- list(...)
  args_forbidden <- args %in% ellipses_forbidden_args
  if(any(args_forbidden)) {
    stop(paste0("args forbidden: ", paste0(unlist(args[args_forbidden]),
                                           collapse = " ")
    ))
  }

  # get the gene universe (all genes in the dataframe)
  gene_universe <- x[[1]]

  # find absolute value of the user-inputted LFC. Accounts for if user inputs a
  # negative LFC.
  lfc_cutoff <- abs(lfc_cutoff)

  # get entrez ids for DE genes, splitting by LFC direction if needed
  if(isTRUE(split_by_lfc_dir)){
    up_genes <- x[(x[[3]] < padj_cutoff) & (x[[2]] > lfc_cutoff),
                  1, drop=TRUE]
    down_genes <- x[(x[[3]] < padj_cutoff) & (x[[2]] < -lfc_cutoff),
                    1, drop=TRUE]
    de_genes <- list(up=up_genes, down=down_genes)

  } else{
    # filter for lfc using abs(LFC)
    de_genes <- x[(x[[3]] < padj_cutoff) &
                    (abs(x[[2]]) > lfc_cutoff), 1, drop=TRUE]

    de_genes <- list(no_dir=de_genes)
  }

  if (identical("reactome", gene_set)){

    enrich_res <- lapply(de_genes, function(gene_list){
      ReactomePA::enrichPathway(gene = gene_list,
                                organism = organism,
                                universe = gene_universe,
                                ...)
    })


  } else if (identical("kegg", gene_set)){
    enrich_res <- lapply(de_genes, function(gene_list){
      clusterProfiler::enrichKEGG(gene = gene_list,
                                  organism = organism,
                                  universe = gene_universe,
                                  ...)
    })

  } else if (gene_set %in% c("H", "C1", "C2", "C3", "C4", "C5", "C6", "C7")){
    msigdbr_cat <- gene_set
    msigdbr_gene_set <- msigdbr::msigdbr(species = organism,
                                         category = gene_set)

    # convert to clusterprofiler friendly format
    msigdbr_clusterprofiler <- msigdbr_gene_set[, c("gs_name", "entrez_gene")] %>%
      as.data.frame(stringsAsFactors = FALSE)

    enrich_res <- lapply(de_genes, function(gene_list){
      clusterProfiler::enricher(gene = gene_list,
                                TERM2GENE = msigdbr_clusterprofiler,
                                universe = gene_universe,
                                ...)
    })

  } else{
    stop("Invalid value for gene_set parameter. See ?run_hypergeometric.")
  }

  enrich_res <- lapply(enrich_res, function(y){
    DOSE::setReadable(enrich_res, OrgDb = orgDb, keyType = "ENTREZID")
  })

  enrich_res
})

###-----------------------------------------------------------------------------
#' @describeIn run_hypergeometric Extract out the gene, LFC and adjusted PValue.
#'   Run run_hypergeometric for each result object (contrast).
#' @export
setMethod("run_hypergeometric", "BbcSE", function(x, de_pkg="edger",
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

    enrich_results <- lapply(edger_results_names, function(edger_res_name){
      dge_table <- edgeR::topTags(edger_results[[edger_res_name]],
                                  n = nrow(dgelist(edger(x)))) %>%
        as.data.frame(stringsAsFactors = FALSE)
      dge_table$entrez_ids <- rowData(x)$entrez_ids[rownames(dge_table)]

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

      # get entrezgene, LFC and adjusted PValue
      genes_lfc_and_adjPval <- filt_dge_table[, c("entrez_ids", "logFC", "FDR")]

      # run the data.frame method of 'run_gsea'
      run_hypergeometric(x = genes_lfc_and_adjPval, ...)
    })

    names(enrich_results) <- edger_results_names

  }

  return(enrich_results)
})
