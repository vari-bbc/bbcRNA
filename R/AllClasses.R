###-----------------------------------------------------------------------------
#' BbcEdgeR, a class for storing an edgeR analysis
#'
#' @importClassesFrom edgeR DGEList DGEGLM DGEExact DGELRT
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom methods slotNames
.BbcEdgeR <- setClass("BbcEdgeR",
                      slots = representation(
                        dgelist = "DGEList",
                        de_results = "list",
                        norm_cts = "SummarizedExperiment"),
                      prototype = prototype(
                        dgelist = edgeR::DGEList(),
                        de_results = list(),
                        norm_cts = SummarizedExperiment()))

setValidity("BbcEdgeR", function(object) {

  msg <- NULL

  # check the slot names
  if (!identical(slotNames(object), c("dgelist","de_results","norm_cts"))) {
    msg <- c(msg, "slot names must be 'dgelist','de_results','norm_cts'")
  }

  # check the dgelist slot
  mydgelist <- dgelist(object)

  if (!is(mydgelist, "DGEList")) {
    msg <- c(msg, "object@dgelist must be a DGEList")
  }

  # check the results slot
  myderesults <- de_results(object)

  if (!is.list(myderesults)) {
    msg <- c(msg, "object@de_results must be a list")
  }

  if(length(myderesults) > 0){
    myderesults_names <- names(myderesults)
    if (!identical(length(myderesults_names),
                   length(unique(myderesults_names)))) {
      msg <- c(msg, "object@de_results must have unique names")
    }

    if (!is(myderesults[[1]], "DGEGLM")) {
      msg <- c(msg,
               "object@de_results[[1]] must be a DGEGLM object.")
    }

    if (!identical(myderesults_names[1], "DGEGLM")) {
      msg <- c(msg,
               "object@de_results[[1]] must be named DGEGLM.")
    }

    invisible(lapply(myderesults[[-1]], function(curr_edger){
      if (!is(curr_edger, "DGEExact") && !is(curr_edger, "DGELRT")) {
        msg <- c(msg,
                 "object@de_results[[-1]] must contain only edgeR result objects.")
      }
    }))

    # check that info in DGEGLM stored in @de_results is same as in @dgelist
    invisible(lapply(c("counts", "samples", "design", "AveLogCPM"), function(x){
      if(!identical(mydgelist[[x]], myderesults$DGEGLM[[x]])){
        msg <- c(msg, "counts, samples, design and AveLogCPM must be the same in dgelist slot and de_results$DGEGLM")
      }
    }))

    # check that info in the edgeR objects in @de_results[[-1]] correspond to that
    # in @de_results[[1]]
    common_elems <- intersect(names(myderesults$DGEGLM), names(myderesults[[2]]))

        if(length(common_elems) != 16) {
      stop("DGEGLM and DGELRT objects in @de_results should share 16 elements with the same name based on test with tcell data")
    }

    for (i in 2:length(myderesults)){
      invisible(lapply(common_elems, function(x){
        if(!identical(myderesults$DGEGLM[[x]], myderesults[[i]][[x]])){
          msg <- c(msg, "Different info in DGEGLM and DGELRT objects in @de_results")
        }
      }))
    }
  }

  # check the norm_cts slot
  mynorm_cts <- norm_cts(object)

  if (!is(mynorm_cts, "SummarizedExperiment")) {
    msg <- c(msg, "object@norm_cts must be a SummarizedExperiment")
  }

  if(nrow(mynorm_cts) > 0){
    # check same genes
    dgelist_genes <- rownames(mydgelist$counts)
    norm_cts_genes <- rownames(mynorm_cts)

    if(!identical(sort(dgelist_genes), sort(norm_cts_genes)) ||
       !identical(length(dgelist_genes), length(norm_cts_genes))){
      msg <- c(msg, "object@norm_cts and object@dgelist must have same genes")
    }

    # check same samples
    dgelist_samples <- colnames(mydgelist$counts)
    norm_cts_samples <- colnames(mynorm_cts)

    if(!identical(sort(dgelist_samples), sort(norm_cts_samples)) ||
       !identical(length(dgelist_samples), length(norm_cts_samples))){
      msg <- c(msg, "object@norm_cts and object@dgelist must have same samples")
    }

    if (assayNames(mynorm_cts)[1] != "norm_log_cpm") {
      msg <- c(msg, "'norm_log_cpm' must be first assay")
    }
  }

  if (is.null(msg)) {
    TRUE
  } else msg
})

###-----------------------------------------------------------------------------
#' BbcSE, an extension of SummarizedExperiment for VAI BBC RNA-seq workflow
#'
#' In a \code{BbcSE} object, "counts" must be the first assay and must contain
#' non-negative values.
#' @importFrom S4Vectors metadata SimpleList
#' @importFrom  SummarizedExperiment assayNames assay
#' @importClassesFrom SummarizedExperiment SummarizedExperiment Assays
#' @importClassesFrom DESeq2 DESeqDataSet
#' @importClassesFrom edgeR DGEList DGEGLM DGEExact DGELRT
.BbcSE <- setClass("BbcSE",
                   slots = representation(
                     aln_rates = "matrix",
                     edger = "BbcEdgeR",
                     deseq2 = "list"),
                   prototype = prototype(
                     assays = Assays(SimpleList(counts=matrix(0, 0, 0))),
                     aln_rates = matrix(0, 0, 0),
                     edger = .BbcEdgeR(),
                     deseq2 = list()),
                   contains="RangedSummarizedExperiment")

setValidity("BbcSE", function(object) {
  NR <- nrow(object)
  NC <- ncol(object)
  msg <- NULL

  if (!identical(length(colnames(object)),
                 length(unique(colnames(object))))) {
    msg <- c(msg, "colnames must be unique")
  }

  if (!identical(length(rownames(object)),
                 length(unique(rownames(object))))) {
    msg <- c(msg, "rownames must be unique")
  }

  if (assayNames(object)[1] != "counts") {
    msg <- c(msg, "'counts' must be first assay")
  }

  if (length(assay(object, "counts")) > 0 && min(assay(object, "counts")) < 0) {
    msg <- c(msg, "negative values in 'counts'")
  }

  ###edger slot-----------------------------------------------------------------
  edger_slot <- edger(object)
  dgelist <- dgelist(edger_slot)

  if (length(dgelist$counts) > 0){
    # check same genes
    if(!all(rownames(dgelist) %in% rownames(object))){
      msg <- c(msg, "Genes in edger slot must be a subset of those in the BbcSE object")
    }

    # check same samples
    if(!identical(sort(colnames(object)), sort(colnames(dgelist)))){
      msg <- c(msg, "Samples must be the same between edger slot and the BbcSE object")
    }

  }


  ###END edger slot-------------------------------------------------------------


  ###deseq2 slot----------------------------------------------------------------
  if (class(deseq2(object))[1] != "list") {
    msg <- c(msg, "deseq2 slot must be a list")
  } else{

    if (length(deseq2(object)) > 0 &&
        !is(deseq2(object)[[1]], "DESeqDataSet")) {
      msg <- c(msg, "deseq2(object)[[1]] must be a DESeqDataSet object")
    }

    if (length(deseq2(object)) > 1) {
      invisible(lapply(deseq2(object)[-1], function(curr_deseq2){
        if (!is(curr_deseq2, "DESeqResults")) {
          msg <- c(msg,
                   "deseq2(object)[-1] elements must be a DESeqResults object.")
        }
      }))
    }
  }
  ###END deseq2 slot------------------------------------------------------------



  ###aln_rates slot-------------------------------------------------------------
  aln_rates <- aln_rates(object, withDimnames=FALSE)

  if (!is.matrix(aln_rates)) {
    msg <- c(msg, "aln_rates must be a matrix")
  } else if (length(aln_rates) > 0){

    aln_rates <- aln_rates(object, withDimnames=TRUE)

    aln_rates_colnames <- colnames(aln_rates)
    if(!identical(length(aln_rates_colnames),
                  length(unique(aln_rates_colnames)))){
      msg <- c(msg, "aln_rates column names must be unique")
    }

    valid_aln_rates_colnames <- c("input_reads",
                                  "uniq_aln_reads",
                                  "mult_aln_reads",
                                  "map_rate",
                                  "uniq_map_rate")
    if(!all(colnames(aln_rates) %in%
           valid_aln_rates_colnames)) {
      msg <- c(
        msg, paste0("colnames for aln_rates must be one of: ",
                    paste(valid_aln_rates_colnames, collapse = ", "))
      )
    }

    if (length(aln_rates(object, withDimnames=FALSE) > 0)) {
      if (nrow(aln_rates(object, withDimnames=FALSE)) != NC) {
        msg <- c(
          msg, "'nrow(aln_rates)' should be equal to the number of columns"
        )
      }

    }
  }
  ###END aln_rates slot---------------------------------------------------------

  if (is.null(msg)) {
    TRUE
  } else msg
})


