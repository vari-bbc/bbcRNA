# test the star_to_mat function ------------------------------------------------
counts_mat <- star_to_mat(dir = ext_data_dir,
                          rgx = "^[^_]+_[^_]+", column = 2)

test_that("star_to_mat returns correct number of columns", {
  expect_equal(length(count_files), ncol(counts_mat))
})

test_that("star_to_mat returns correct number of rows", {
  column_names <- c("gene_id", "unstr", "str_r1", "str_r2")

  # read in first file
  counts_df <- readr::read_tsv(file = count_files[1],
                               col_names = column_names, col_types = "ciii")

  # In typical STAR output, 4 rows contain non-gene data. They begin with "N_"
  counts_df <- dplyr::filter(counts_df, !stringr::str_detect(gene_id, "^N_"))

  # Now, number of rows is the number of genes in the first file
  ngenes_first_cts_file <- nrow(counts_df)

  expect_equal(ngenes_first_cts_file, nrow(counts_mat))
})

# test the read_star_map_rates function ----------------------------------------
aln_rates <- read_star_map_rates(dir = ext_data_dir,
                                 rgx = "^[^_]+_[^_]+")

test_that("read_star_map_rates returns correct number of rows", {
  expect_equal(nrow(aln_rates), length(aln_metrics_files))
})

# test the read_col_meta function ----------------------------------------------
col_meta <- read_col_meta(paste0(ext_data_dir, "/meta.txt"))

test_that("read_col_meta returns DataFrame", {
  expect_s4_class(col_meta, "DataFrame")
})

test_that("read_col_meta returns DataFrame with rownames", {
  expect_gt(length(rownames(col_meta)), 0)
})



