# test the get_de_results function ---------------------------------------------
test_that("get_de_results returns correct data", {
  # test glmQLFTest output
  out_data_qlf <- get_de_table(bbc_obj_glmQLFTest, de_pkg="edger")[[1]]

  temp_edger <- edger(bbc_obj_glmQLFTest)

  ## check that all genes are returned
  expect_equal(nrow(out_data_qlf), nrow(dgelist(temp_edger)))

  de_results_qlf <- as.data.frame(edgeR::topTags(
    de_results(edger(bbc_obj_glmQLFTest))[[2]],
    n=nrow(dgelist(edger(bbc_obj_glmQLFTest)))))

  de_results_qlf$gene_symbol <-
    rowData(bbc_obj_glmQLFTest)[rownames(de_results_qlf), "uniq_syms"]

  group_mean_expr <- rowData(norm_cts(temp_edger)) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::select(dplyr::ends_with("norm_log_cpm"))

  de_results_qlf <-
    cbind(de_results_qlf,
          group_mean_expr[rownames(de_results_qlf), ])

  de_results_qlf <- de_results_qlf %>%
    tibble::rownames_to_column("ensembl_id")

  ## check that same DE results returned
  expect_equal(ncol(out_data_qlf), ncol(de_results_qlf))
  expect_equal(out_data_qlf, de_results_qlf[, colnames(out_data_qlf)])

  ### --------------------------------------------------------------------------

  # test glmTreat output
  out_data_glmtreat <- get_de_table(bbc_obj_glmTreat, de_pkg="edger")[[1]]

  temp_edger <- edger(bbc_obj_glmTreat)

  ## check that all genes are returned
  expect_equal(nrow(out_data_glmtreat), nrow(dgelist(temp_edger)))

  de_results_glmtreat <- as.data.frame(edgeR::topTags(
    de_results(edger(bbc_obj_glmTreat))[[2]],
    n=nrow(dgelist(edger(bbc_obj_glmTreat)))))

  de_results_glmtreat$gene_symbol <-
    rowData(bbc_obj_glmTreat)[rownames(de_results_glmtreat), "uniq_syms"]

  group_mean_expr <- rowData(norm_cts(temp_edger)) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    dplyr::select(dplyr::ends_with("norm_log_cpm"))

  de_results_glmtreat <-
    cbind(de_results_glmtreat,
          group_mean_expr[rownames(de_results_glmtreat), ])

  de_results_glmtreat <- de_results_glmtreat %>%
    tibble::rownames_to_column("ensembl_id")

  ## check that same DE results returned
  expect_equal(ncol(out_data_glmtreat), ncol(de_results_glmtreat))
  expect_equal(out_data_glmtreat,
               de_results_glmtreat[, colnames(out_data_glmtreat)])

})
