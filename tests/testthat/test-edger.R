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

  temp_dgelist <- temp_dgelist[rows_keep, , keep.lib.sizes=FALSE]

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

  expect_equivalent(qlf_test, de_results(edger(bbc_obj_glmQLFTest))[[2]])
  expect_equivalent(glmtreat_test, de_results(edger(bbc_obj_glmTreat))[[2]])
})

test_that("findDEGs does not produce an error with two contrasts", {
  # should not produce an error
  expect_error(findDEGs(bbc_obj_dgelist,
                        de_pkg = "edger",
                        test = "glmQLFTest",
                        design = "~0+genotype",
                        contrasts = list(c("genotype", "mut", "WT"),
                                         c("genotype", "WT", "mut"))), NA)
})

test_that("findDEGs with coefs returns same table as contrasts", {
  expect_equivalent(de_results(edger(bbc_obj_glmQLFTest_coef))[[2]]$table,
                    de_results(edger(bbc_obj_glmQLFTest))[[2]]$table)
  expect_equivalent(de_results(edger(bbc_obj_glmTreat_coef))[[2]]$table,
                    de_results(edger(bbc_obj_glmTreat))[[2]]$table)
})

test_that("findDEGs works when both coefs and contrasts supplied", {

  combined <- findDEGs(x=bbc_obj_dgelist, design="~0+genotype",
                       coefs=list(1, 2),
                       contrasts=list(c("genotype", "mut", "WT")))

  expect_equivalent(de_results(edger(combined))[["genotype_mut_-_WT"]],
                    de_results(edger(bbc_obj_glmQLFTest))[[2]])

  coef <- findDEGs(x=bbc_obj_dgelist, design="~0+genotype",
                       coefs=list(2))

  expect_equivalent(de_results(edger(combined))[["coef_2"]],
                    de_results(edger(coef))[[2]])
})

test_that("findDEGs works when both coefs specified as name", {

  coef <- findDEGs(x=bbc_obj_dgelist, design="~genotype",
                   coefs=list("genotypemut"))

  expect_equivalent(de_results(edger(bbc_obj_glmQLFTest_coef))[[2]],
                    de_results(edger(coef))[[2]])
})

test_that("findDEGs works when multiple coefs specified", {

  coef <- findDEGs(x=bbc_obj_dgelist, design="~genotype",
                   coefs=list(2, 1:2))

  expect_equivalent(de_results(edger(bbc_obj_glmQLFTest_coef))[[2]],
                    de_results(edger(coef))[[2]])

  design <- model.matrix(as.formula("~genotype"),
                         data = colData(bbc_obj_dgelist))

  temp_dgelist <- dgelist(edger(bbc_obj_dgelist))

  # calculate dispersions
  temp_dgelist <- edgeR::estimateDisp(temp_dgelist, design, robust = TRUE)

  # calculate fit
  temp_fit <- edgeR::glmQLFit(temp_dgelist, design, robust = TRUE)

  qlf_test <- edgeR::glmQLFTest(temp_fit,
                                coef = 1:2)

  expect_equivalent(qlf_test,
                    de_results(edger(coef))[[3]])
})

test_that("Invalid coefs return error", {
  expect_error(findDEGs(x=bbc_obj_dgelist, design="~genotype",
                        coefs=list("foobar")),
               regexp="Invalid coefficient: foobar")
  expect_error(findDEGs(x=bbc_obj_dgelist, design="~genotype",
                        coefs=list(100)),
               regexp="Invalid coefficient: 100")
})


