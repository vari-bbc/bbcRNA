#' Convert STAR read counts to count matrix.
#'
#' \code{star_to_mat} returns a count matrix from STAR ReadsPerGene.out.tab
#' files.
#'
#' This function will read in all the ReadsPerGene.out.tab files in the
#' specified directory and converts it to a counts matrix.
#'
#' @param dir A character scalar indicating the directory containing all the
#'     STAR ReadsPerGene.out.tab files.
#' @param rgx A character scalar representing a regex used to parse out the
#'     sample name from the name of the ReadsPerGene.out.tab file.
#' @param column An integer indicating the column to extract counts from.
#'     1 = unstranded, 2 = 1st read strand aligned with RNA, 3 = 2nd read
#'     strand aligned with RNA  with RNA
#' @param rm_ens_vers Logical indicating whether version number should be
#'     removed from gene ID (as indicated by trailing period and integers).
#' @return A matrix where the column names are the sample names and the row
#'     names are the gene names
#' @export
star_to_mat <- function(dir, rgx, column, rm_ens_vers = TRUE){

  # get all the file names for the count files from STAR
  count_files <- grep("*ReadsPerGene.out.tab",
                      list.files(dir, full.names = TRUE),
                      value = TRUE)

  stopifnot(length(count_files) > 0) # check at least 1 count file

  stopifnot(column %in% c(1, 2, 3)) # make sure column requested is valid

  # read in files
  counts_list <- lapply(count_files,
                        function(x){
                          column_names <- c("gene_id",
                                            "unstr",
                                            "str_r1",
                                            "str_r2")

                          # read
                          df <- readr::read_tsv(file = x,
                                                col_names = column_names,
                                                col_types = "ciii")

                          # add 'sample' column based on filename
                          samp_name <- stringr::str_extract(basename(x), rgx)

                          df$sample = samp_name

                          return(df)
                        })

  # check same columns in every file
  colnames_list <- lapply(counts_list, colnames)
  stopifnot(all(sapply(colnames_list[-1], identical, colnames_list[[1]])))

  # get the counts column indicated by function call
  counts_col_name <- colnames_list[[1]][1 + column]

  # process files
  counts_list <- lapply(counts_list, function(x){

    df <- x %>%
      dplyr::filter(!.data$gene_id %in% c("N_unmapped",
                                          "N_multimapping",
                                          "N_noFeature",
                                          "N_ambiguous")) %>%
      dplyr::select(.data$gene_id,
                    dplyr::matches(counts_col_name),
                    .data$sample)

    # check for missing data
    stopifnot(all(complete.cases(df)))

    return(df)

  })

  # check each file has same # of genes
  samp_lens <- unlist(lapply(counts_list, function(x){
    num_row <- nrow(x)
    num_uniq_genes <- length(unique(x$gene_id))
    stopifnot(num_row == num_uniq_genes) # check gene names are unique

    return(num_row)
    }))
  samp_lens_uniq <- unique(unname(samp_lens))
  stopifnot(length(samp_lens_uniq) == 1 & samp_lens_uniq > 1)

  # rbind
  counts_rbind <- do.call(rbind, counts_list)

  # check each file has same genes
  ## split by gene
  split_by_genes <- split(counts_rbind, counts_rbind$gene_id)

  ## samples per gene
  split_by_genes_nrows <- unlist(lapply(split_by_genes, nrow))

  ## unique # of samples per gene
  split_by_genes_nrows_uniq <- unique(unname(split_by_genes_nrows))

  ## check that all genes have same number of samples and that # is equal
  ## to the number of input samples.
  stopifnot(length(split_by_genes_nrows_uniq) == 1 &
              split_by_genes_nrows_uniq == length(count_files))

  # spread to form count matrix
  count_mat <- tidyr::spread(data = counts_rbind,
                             key = .data$sample,
                             value = !!counts_col_name) %>%
    tibble::column_to_rownames(var = "gene_id")

  # check # of genes unchanged
  stopifnot(nrow(count_mat) == samp_lens_uniq)

  # remove version number from gene ids if rm_ens_vers
  if (rm_ens_vers) {
    new_gene_ids <- stringr::str_remove(rownames(count_mat), "\\.\\d+$")

    # check that removing the version # does not affect uniqueness.
    stopifnot(length(unique(new_gene_ids)) == length(rownames(count_mat)))

    # set new rownames
    rownames(count_mat) <- new_gene_ids
  }

  count_mat <- as.matrix(count_mat)

  return(count_mat)
}


#' Import STAR Log.final.out files for analyzing mapping rates.
#'
#' \code{read_star_map_rates} returns a dataframe from STAR Log.final.out files.
#'
#' This function will read in all the Log.final.out files in the specified
#' directory and returns a dataframe with alignment metrics.
#'
#' @param dir A character scalar indicating the directory containing all the
#'     STAR Log.final.out files.
#' @param rgx A character scalar representing a regex used to parse out the
#'     sample name from the name of the Log.final.out files.
#' @return A dataframe.
#' @export
read_star_map_rates <- function(dir, rgx){

  # get all the file names for the aln metrics files from STAR
  aln_metrics_files <- grep(".*Log.final.out$",
                      list.files(dir, full.names = TRUE),
                      value = TRUE,
                      perl = TRUE)

  stopifnot(length(aln_metrics_files) > 0) # check at least 1 file

  # read in files
  metrics_of_interest <- c(input_reads =
                             "Number of input reads",
                           uniq_aln_reads =
                             "Uniquely mapped reads number",
                           mult_aln_reads =
                             "Number of reads mapped to multiple loci")

  # parse out rows of interest
  aln_metrics_list <- lapply(aln_metrics_files,
                             function(x){
                               samp_name <- stringr::str_extract(basename(x),
                                                                 rgx)

                               metrics <- readr::read_delim(x,
                                                 delim = "\\n",
                                                 col_names = FALSE,
                                                 col_types = "c")

                               stopifnot(ncol(metrics) == 1)

                               keep_rows_list <- lapply(metrics_of_interest,
                                                        grep,
                                                        x = metrics[[1]])

                               keep_rows <- metrics[unlist(keep_rows_list), ]

                               keep_metrics <- stringr::str_extract(
                                 keep_rows[[1]],
                                 "(?<=\\\t)\\d+$")

                               names(keep_metrics) <- names(metrics_of_interest)

                               keep_metrics["sample"] <- samp_name

                               keep_metrics_df <- data.frame(
                                 as.list(keep_metrics),
                                 stringsAsFactors = FALSE)

                               return(keep_metrics_df)
                             })

  aln_metrics_df <- do.call(rbind, aln_metrics_list)

  return(aln_metrics_df)
}
