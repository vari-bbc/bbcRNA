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

  return(count_mat)
}

