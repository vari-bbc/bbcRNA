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
#'   removed. Implements edgeR::filterByExpr. Filters DGEList only.
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
  if(nrow(dgelist(edger(x))) > 0) stop("edger slot not empty")
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
  edger(x) <- BbcEdgeR(dgelist = myDGEList)

  # add log normalized counts if norm factors calculated
  if (isTRUE(calc_norm)) {
    x <- normalize_counts(x)
  }

  validObject(x)

  return(x)
}

###-----------------------------------------------------------------------------
#' @describeIn normalize_counts For de_pkg="edger", a SummarizedExperiment
#'   object is created and stored in the BbcEdgeR object in the edger slot.
#'   Group average log normalized counts are also calculated and stored as
#'   rowData.
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment assay
#' @export
setMethod("normalize_counts", "BbcSE", function(x, de_pkg = "edger") {

  if(!de_pkg %in% c("edger", "deseq2")){
    stop("de_pkg not recognized")
  }

  if(identical(de_pkg, "edger")){
    dgelist_obj <- dgelist(edger(x))
    if(nrow(dgelist_obj) == 0){
      stop("DGEList in edger slot is empty. Run makeDGEList first.")
    }

    # calculate normalized counts. Use normalize_counts,DGEList
    norm_log_cpm <- normalize_counts(dgelist_obj)

    # calculate by group cpms
    group_norm_log_cpm <- edgeR::cpmByGroup(dgelist_obj,
                                            normalized.lib.sizes = TRUE,
                                            log=TRUE)

    ## modify col names of group cpms
    colnames(group_norm_log_cpm) <- paste0(colnames(group_norm_log_cpm),
                                           ".norm_log_cpm")

    ## add group cpms to rowData.
    new_row_data <- cbind(rowData(x)[rownames(dgelist_obj), , drop = FALSE],
                          group_norm_log_cpm)

    # store normalized counts in a SummarizedExperiment
    se <- SummarizedExperiment(
      assays = list(norm_log_cpm = norm_log_cpm),
      rowData = new_row_data,
      colData = colData(x)
    )

  }

  # extract the BbcEdgeR object from the edger slot of the BbcSE object
  bbcedger_obj <- edger(x)

  # store the normalized counts
  norm_cts(bbcedger_obj) <- se

  # replace the BbcEdgeR object in the edger slot of the BbcSE object
  edger(x) <- bbcedger_obj

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

###-----------------------------------------------------------------------------

#' Function to convert to edgeR-style contrasts
#'
#' @param contrasts_list list of chr vectors containing variable name, numerator
#'     level, denominator level
#' @param design output of model.matrix
#' @importFrom stringr str_detect str_replace
#' @importFrom limma makeContrasts
#' @export
lists_to_edger_contrasts <- function(contrasts_list, design){
  out_names <- sapply(contrasts_list, function(x){
    stopifnot(length(x) == 3 && is.vector(x, mode="character"))
    paste0(x[1], "_", x[2], "_-_", x[3])
  })
  out_contrasts <- lapply(contrasts_list, function(x){
    stopifnot(length(x) == 3 && is.vector(x, mode="character"))

    de_var <- x[1]

    num_and_denom <- x[2:3]

    num_and_denom <- sapply(num_and_denom, function(y){

      # splits if numerator was a difference between two levels
      # (ie. interaction)
      y <- stringr::str_split(y, "-", simplify = TRUE)[1,]

      prefixed <- sapply(y, function(z) paste0(de_var, z))

      if(length(prefixed) > 1){
        prefixed <- paste0(prefixed[1], "-", prefixed[2])
      }

      prefixed
    })

    paste0("(", num_and_denom[1], ")-(", num_and_denom[2], ")")
  })

  edger_contrasts <- do.call(function(...){
    limma::makeContrasts(..., levels = design)
  }, out_contrasts)

  colnames(edger_contrasts) <- out_names

  return(edger_contrasts)
}


###-----------------------------------------------------------------------------
#' @describeIn findDEGs Run glmQLFTest or glmTreat. Accepts contrasts in the
#'   format (a matrix) produced by limma::makeContrasts
#' @importFrom edgeR glmQLFTest glmTreat glmQLFit estimateDisp
#' @importFrom limma nonEstimable
#' @importFrom stats model.matrix as.formula
#' @importFrom stringr str_detect
#' @export
setMethod("findDEGs", "DGEList", function(x, test, design, contrasts, coefs,
                                          sample_meta, lfc = log2(2)) {

  design <- model.matrix(as.formula(design), data = sample_meta)

  # dashes will be confusing when defining 'interaction' contrasts.
  if (any(stringr::str_detect(colnames(design), "-"))) {
    stop("No dashes allowed in design colnames")
  }

  # check if anything in the design is non-estimable:
  if (!is.null(limma::nonEstimable(design))) stop("Design not estimable.")

  # calculate dispersions
  x <- edgeR::estimateDisp(x, design, robust = TRUE)

  # calculate fit
  fit <- edgeR::glmQLFit(x, design, robust = TRUE)

  test_res_contrasts <- list()
  test_res_coefs <- list()

  if (!is.null(contrasts)){
    # Format contrasts so that they are readable for edgeR
    edger_contrasts <- lists_to_edger_contrasts(contrasts_list = contrasts,
                                                design = design)
    # either test produces a DGELRT object
    if (identical(test, "glmQLFTest")){

      test_res_contrasts <- lapply(
        colnames(edger_contrasts), function(contr_name) {
          edgeR::glmQLFTest(fit, contrast = edger_contrasts[, contr_name])
        }
      )

    } else if (identical(test, "glmTreat")){

      test_res_contrasts <- lapply(
        colnames(edger_contrasts), function(contr_name) {
          edgeR::glmTreat(fit, contrast = edger_contrasts[, contr_name],
                          lfc = lfc)
        }
      )
    }
    # set the contrasts as the names for each contrast in test_res_contrasts
    names(test_res_contrasts) <- colnames(edger_contrasts)
  }

  if (!is.null(coefs)){
    # check coefs valid
    for (coef in coefs){
      # use sapply to account for when testing multiple coefs
      tryCatch({
        stopifnot(all(sapply(as.data.frame(design)[coef], is.numeric)))
      }, error = function(e) {
        stop(paste0("Invalid coefficient: ", coef))
      })
    }

    # either test produces a DGELRT object
    if (identical(test, "glmQLFTest")){

      test_res_coefs <- lapply(
        coefs, function(coef) {
          edgeR::glmQLFTest(fit, coef = coef)
        }
      )

    } else if (identical(test, "glmTreat")){

      test_res_coefs <- lapply(
        coefs, function(coef) {
          edgeR::glmTreat(fit, coef = coef, lfc = lfc)
        }
      )

    }
    # set the coefs as the names for each coef in test_res_coefs
    # use make.names in case there is a ":" (testing multiple coefs at the same time)
    names(test_res_coefs) <- make.names(paste0("coef_", unlist(coefs)))
  }


  # store the fit in the first element of the results list and the results from
  # each contrast in the rest of the list
  edger_res_list <- c(list(DGEGLM = fit), test_res_contrasts, test_res_coefs)

  return(BbcEdgeR(dgelist = x, de_results = edger_res_list))
})

###-----------------------------------------------------------------------------
#' @describeIn findDEGs Run edgeR or DESeq2 DE testing workflows.
#' @importFrom edgeR glmQLFTest glmTreat glmQLFit
#' @export
setMethod("findDEGs", "BbcSE", function(x, de_pkg = "edger",
                                        test = "glmQLFTest", design,
                                        contrasts = NULL,
                                        coefs = NULL,
                                        lfc = log2(2)) {

  if((!is.null(contrasts) & !all(!duplicated(contrasts))) |
     (!is.null(coefs) & !all(!duplicated(coefs)))){
    stop("Contrasts and coefs cannot be duplicated.")
  }

  if(is.null(contrasts) & is.null(coefs)){
    stop("Specify contrasts and/or coefs")
  }

  if (identical(de_pkg, "edger")){
    if("results" %in% names(edger(x))) {
      stop("edger slot already contains results")
    }

    bbcedger_obj <- findDEGs(x = dgelist(edger(x)), test = test,
                             design = design, contrasts = contrasts,
                             coefs = coefs, sample_meta = colData(x), lfc = lfc)

    # copy norm_cts slot from the old BbcEdgeR object
    norm_cts(bbcedger_obj) <- norm_cts(edger(x))

    # replace the edger slot with the new BbcEdgeR object containing the results
    edger(x) <- bbcedger_obj
  }

  validObject(x)

  return(x)
})

