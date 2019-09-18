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
#' @return A matrix where the column names are the sample names and the row
#'     names are the gene names
#' @export
star_to_mat <- function(dir, rgx, column){

  # get all the file names for the count files from STAR
  count_files <- grep("*ReadsPerGene.out.tab",
                      list.files(dir),
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
                                                col_names = column_names)

                          # add 'sample' column based on filename
                          samp_name <- stringr::str_extract(x, rgx)

                          df$sample = samp_name

                          return(df)
                        })

  # process files
  counts_list <- lapply(counts_list, function(x){
    # get the counts column indicated by function call
    counts_col_name <- colnames(x)[1+column]

    x %>%
      dplyr::filter(!.data$gene_id %in% c("N_unmapped",
                                          "N_multimapping",
                                          "N_noFeature",
                                          "N_ambiguous")) %>%
      dplyr::select(.data$gene_id,
                    dplyr::matches(counts_col_name),
                    .data$sample)

  })

  # rbind
  counts_rbind <- do.call(rbind, counts_list)

  # spread to form count matrix
  # note that 'spread_' means standard evaluation
  count_mat <- counts_rbind %>%
    tidyr::spread(.data, key = 3, value = 2) %>%
    tibble::column_to_rownames(.data, var = "gene_id")

  return(count_mat)
}

