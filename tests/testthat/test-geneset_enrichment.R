test_that("run_gsea runs without errors for various ranking metrics", {
  expect_error(run_gsea(bbc_obj_glmQLFTest, gene_set="H",
                        organism="Mus musculus", rank_by="signed-log10pval",
                        orgDb=org.Mm.eg.db), NA)
  expect_error(run_gsea(bbc_obj_glmQLFTest, gene_set="H",
                        organism="Mus musculus", rank_by="-log10pval",
                        orgDb=org.Mm.eg.db), NA)
})
