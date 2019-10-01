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
