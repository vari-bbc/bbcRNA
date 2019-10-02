# test the makeDGEList function ------------------------------------------------
test_that("makeDGEList runs correctly", {
  # make DGEList object, calculate normalization factors
  bbc_obj2a <- makeDGEList(bbc_obj, group="genotype")

  # Run edgeR methods outside of function
  # make a copy of bbc_obj
  bbc_obj3 <- bbc_obj

  # make DGEList
  myDGEList <- edgeR::DGEList(counts = assay(bbc_obj3, "counts"),
                              group = colData(bbc_obj3)[["genotype"]])

  # filter out lowly expressed genes
  rows_keep <- edgeR::filterByExpr(myDGEList,
                                   group = myDGEList$samples$group)

  myDGEList <- myDGEList[rows_keep, ]

  # calculate norm factors
  myDGEList <- edgeR::calcNormFactors(myDGEList)

  # store DGEList in the BbcSE object
  edger(bbc_obj3) <- BbcEdgeR(dgelist = myDGEList)

  bbc_obj3 <- normalize_counts(bbc_obj3)

  validObject(bbc_obj3)

  expect_equivalent(edger(bbc_obj2a), edger(bbc_obj3))

  # test rm_low_genes = FALSE -----------------------------------------------------

  bbc_obj2b <- makeDGEList(bbc_obj, group="genotype", rm_low_genes=FALSE)
  bbc_obj2b_norm_cts <- norm_cts(edger(bbc_obj2b))

  # remove normalization factors and compare with makeDGEList(bbc_obj,
  # group="genotype", calc_norm = FALSE)
  temp_dgelist <- edgeR::DGEList(counts = assay(bbc_obj3, "counts"),
                                 group = colData(bbc_obj3)[["genotype"]])

  # calculate norm factors
  temp_dgelist <- edgeR::calcNormFactors(temp_dgelist)

  # normalized counts
  norm_cts <- (normalize_counts(temp_dgelist))

  # calculate by group cpms
  group_norm_log_cpm <- edgeR::cpmByGroup(temp_dgelist,
                                          normalized.lib.sizes = TRUE,
                                          log=TRUE)

  ## modify col names of group cpms
  colnames(group_norm_log_cpm) <- paste0(colnames(group_norm_log_cpm),
                                         ".norm_log_cpm")

  ## add group cpms to rowData.
  new_row_data <- cbind(rowData(bbc_obj3)[rownames(temp_dgelist), , drop = FALSE],
                        group_norm_log_cpm)

  # store normalized counts in a SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(norm_log_cpm = norm_cts),
    rowData = new_row_data,
    colData = colData(bbc_obj3)
  )

  expect_equivalent(dgelist(edger(bbc_obj2b)), temp_dgelist)
  expect_equivalent(assay(bbc_obj2b_norm_cts), norm_cts)
  expect_equivalent(rowData(bbc_obj2b_norm_cts), new_row_data)
  expect_equivalent(colData(bbc_obj2b_norm_cts), colData(bbc_obj3))
  expect_equivalent(bbc_obj2b_norm_cts, se)


  # test calc_norm = FALSE -----------------------------------------------------

  bbc_obj2c <- makeDGEList(bbc_obj, group="genotype", calc_norm = FALSE)

  # remove normalization factors and compare with makeDGEList(bbc_obj,
  # group="genotype", calc_norm = FALSE)
  temp_dgelist <- edgeR::DGEList(counts = assay(bbc_obj3, "counts"),
                                 group = colData(bbc_obj3)[["genotype"]])
  temp_dgelist$samples$norm.factors <- 1

  # filter out lowly expressed genes
  rows_keep <- edgeR::filterByExpr(temp_dgelist,
                                   group = temp_dgelist$samples$group)

  temp_dgelist <- temp_dgelist[rows_keep, ]

  expect_equivalent(dgelist(edger(bbc_obj2c)), temp_dgelist)
  expect_equivalent(norm_cts(edger(bbc_obj2c)), SummarizedExperiment())
})

# test the findDEGs function ---------------------------------------------------
test_that("findDEGs runs correctly", {
  design <- model.matrix(as.formula("~0+genotype"),
                         data = colData(bbc_obj_dgelist))

  contrast <- limma::makeContrasts("genotypemut-genotypeWT", levels = design)

  temp_dgelist <- dgelist(edger(bbc_obj_dgelist))

  # calculate dispersions
  temp_dgelist <- edgeR::estimateDisp(temp_dgelist, design, robust = TRUE)

  # calculate fit
  temp_fit <- edgeR::glmQLFit(temp_dgelist, design, robust = TRUE)

  qlf_test <- edgeR::glmQLFTest(temp_fit,
                                contrast = contrast)

  glmtreat_test <- edgeR::glmTreat(temp_fit,
                                   contrast = contrast,
                                   lfc = log2(2))

  ## Identify DE genes with glmQLFTest
  bbc_obj_glmQLFTest <- findDEGs(bbc_obj_dgelist,
                                 de_pkg = "edger",
                                 test = "glmQLFTest",
                                 design = "~0+genotype",
                                 contrasts = list(c("genotype", "mut", "WT")))

  ## Identify DE genes with glmTreat
  bbc_obj_glmTreat <- findDEGs(bbc_obj_dgelist,
                               de_pkg = "edger",
                               test = "glmTreat",
                               design = "~0+genotype",
                               contrasts = list(c("genotype", "mut", "WT")))

  str(qlf_test)
  str(de_results(edger(bbc_obj_glmQLFTest))[[2]])
  expect_equivalent(qlf_test, de_results(edger(bbc_obj_glmQLFTest))[[2]])
  expect_equivalent(glmtreat_test, de_results(edger(bbc_obj_glmTreat))[[2]])
})


