test_that("run_gsea runs without errors for various ranking metrics", {
  expect_error(run_gsea(bbc_obj_glmQLFTest, gene_set="H",
                        organism="Mus musculus", rank_by="signed-log10pval",
                        orgDb=org.Mm.eg.db), NA)
  expect_error(run_gsea(bbc_obj_glmQLFTest, gene_set="H",
                        organism="Mus musculus", rank_by="-log10pval",
                        orgDb=org.Mm.eg.db), NA)
})

test_that("find_missing_in_gseaResult results consistent with run_gsea messages in terms of genes removed", {

  gsea_results_list <- run_gsea(bbc_obj_glmQLFTest, gene_set="H",
                                organism="Mus musculus", rank_by="signed-log10pval",
                                orgDb=org.Mm.eg.db)

  msges <- capture_messages(run_gsea(bbc_obj_glmQLFTest, gene_set="H",
                                     organism="Mus musculus", rank_by="signed-log10pval",
                                     orgDb=org.Mm.eg.db))

  missing_genes <- find_missing_in_gseaResult(gsea_results_list, bbc_obj_glmQLFTest)


  gsea_results_genes_rmed <- lapply(names(gsea_results_list), function(gsea_res_name){
    curr_msges <- grep(gsea_res_name, msges, value=TRUE) # grep for message corresponding to specific contrast
    genes_rm_due_to_lfc_or_pval <- stringr::str_extract(string = curr_msges[1], pattern="(?<=removed )\\d+(?= genes)") %>%
      as.integer()
    genes_rm_due_to_no_entrez_match <- stringr::str_extract(string = curr_msges[2], pattern="(?<=removed )\\d+(?= genes)") %>%
      as.integer()
    expect_equivalent(genes_rm_due_to_lfc_or_pval, missing_genes[[gsea_res_name]] %>%
                        dplyr::filter(logFC==0 | PValue==1) %>%
                        nrow())
    expect_equivalent(genes_rm_due_to_no_entrez_match, missing_genes[[gsea_res_name]] %>%
                        dplyr::filter(low_expr==FALSE & is.na(entrez_ids)) %>%
                        nrow())

  })

})
