# test the plot_aln_metrics function -------------------------------------------
test_that("plot_aln_metrics returns correct data", {
  plot <- plot_aln_metrics(bbc_obj)
  aln_metrics <- aln_metrics(bbc_obj)
  sample_meta <- colData(bbc_obj)
  combined_df <- cbind(aln_metrics, sample_meta[rownames(aln_metrics), ])
  combined_df <- combined_df %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(input_reads = round(.data$input_reads / 10^6, digits=2),
                  uniq_aln_reads = round(.data$uniq_aln_reads / 10^6, digits=2),
                  mult_aln_reads = round(.data$mult_aln_reads / 10^6, digits=2))

  expect_equivalent(plot$data, combined_df[, colnames(plot$data)])
})

# plot_heatmap -----------------------------------------------------------------
test_that("plot_heatmap works for one gene", {
  # expect no error
  expect_error(
    plot_heatmap(bbc_obj_glmQLFTest,
                 genes = rownames(bbc_obj_glmQLFTest@edger@dgelist)[1],
                 zscores = TRUE,
                 gene_labels = "uniq_syms"), NA)
})

test_that("plot_heatmap works when there are requested genes with no normalized counts", {
  genes_in_SE <- nrow(bbc_obj_glmQLFTest)
  genes_in_dgelist <- nrow(bbc_obj_glmQLFTest@edger@dgelist)
  genes_rm <- genes_in_SE - genes_in_dgelist
  msg1 <- paste0("Not plotting ", genes_rm, " genes because normalized counts not available.")
  msg2 <- paste0("Plotting ", genes_in_dgelist, " genes with normalized counts.")

  msges <- capture_messages(plot_heatmap(bbc_obj_glmQLFTest, zscores = FALSE))
  expect_match(msges[1], msg1)
  expect_match(msges[2], msg2)

  expect_error(plot_heatmap(bbc_obj_glmQLFTest, genes = c("foo", "bar")),
  "No genes with normalized counts available.")
})
