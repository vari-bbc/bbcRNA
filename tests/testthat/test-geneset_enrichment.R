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

test_that("run_hypergeometric retruns correct results", {
  contrast_DGELRT <- bbc_obj_glmQLFTest@edger@de_results[[2]] # first element is a DGEGLM
  contrast_DGELRT_toptags <- edgeR::topTags(contrast_DGELRT,
                                            n = nrow(contrast_DGELRT$table)) %>%
    as.data.frame(stringsAsFactors=TRUE)
  all_ens <- rownames(contrast_DGELRT_toptags)
  all_entrez <- rowData(bbc_obj_glmQLFTest)[all_ens, "entrez_ids", drop=TRUE]
  all_entrez <- na.omit(all_entrez)

  sig_up_ens <- rownames(contrast_DGELRT_toptags)[contrast_DGELRT_toptags$FDR < 0.05 & contrast_DGELRT_toptags$logFC > 0]
  sig_up_entrez <- rowData(bbc_obj_glmQLFTest)[sig_up_ens, "entrez_ids", drop=TRUE]
  sig_up_entrez <- na.omit(sig_up_entrez)

  sig_down_ens <- rownames(contrast_DGELRT_toptags)[contrast_DGELRT_toptags$FDR < 0.05 & contrast_DGELRT_toptags$logFC < 0]
  sig_down_entrez <- rowData(bbc_obj_glmQLFTest)[sig_down_ens, "entrez_ids", drop=TRUE]
  sig_down_entrez <- na.omit(sig_down_entrez)

  # the first element of the results of 'run_hypergeometric' here should be the
  # enrichment results for the up-regulated genes for the first contrast
  expect_equivalent(run_hypergeometric(bbc_obj_glmQLFTest, gene_set="kegg",
                                       organism="mmu", orgDb=org.Mm.eg.db)[[1]],
                    clusterProfiler::enrichKEGG(gene=sig_up_entrez,
                                                universe=all_entrez,
                                                organism="mmu"))

  expect_equivalent(run_hypergeometric(bbc_obj_glmQLFTest, gene_set="kegg",
                                       organism="mmu", orgDb=org.Mm.eg.db)[[2]],
                    clusterProfiler::enrichKEGG(gene=sig_down_entrez,
                                                universe=all_entrez,
                                                organism="mmu"))

  expect_equivalent(run_hypergeometric(bbc_obj_glmQLFTest, gene_set="reactome",
                                       organism="mouse", orgDb=org.Mm.eg.db)[[1]],
                    ReactomePA::enrichPathway(gene=sig_up_entrez,
                                              universe=all_entrez,
                                              organism="mouse"))

})
