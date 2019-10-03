#' Get a table of DE results
#'
#' @param x A BbcSE object
#' @param de_pkg "edger" or "deseq2"
#' @return A list of datafrane for each contrast.
#' @importFrom dplyr select ends_with left_join arrange rename
#' @importFrom tibble rownames_to_column
#' @importFrom edgeR topTags
#' @export
get_de_table <- function(x, de_pkg = "edger"){
  if(!is(x, "BbcSE")){
    stop("x must be a BbcSE object.")
  }

  if(identical(de_pkg, "edger")){
    myBbcEdgerobj <- edger(x)

    # get the contents from the de_results slot, removing the fit object (1st
    # element)
    my_de_results <- de_results(myBbcEdgerobj)[-1]

    # get names of all contrasts
    my_edger_res_names <- names(my_de_results)

    res_list <- lapply(my_edger_res_names, function(res_name){
      top_tags <- edgeR::topTags(my_de_results[[res_name]],
                               n = nrow(dgelist(myBbcEdgerobj))) %>%
        as.data.frame(stringsAsFactors = FALSE)

      top_tags_colnames <- colnames(top_tags)

      top_tags <- top_tags %>%
        tibble::rownames_to_column("ensembl_id")

      # get the gene symbols
      if (!"uniq_syms" %in% colnames(rowData(x))){
        stop("Please run 'ens2syms' first.")
      }
      gene_symbols <- rowData(x)[, "uniq_syms", drop = FALSE] %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        tibble::rownames_to_column("ensembl_id") %>%
        dplyr::rename("gene_symbol" = "uniq_syms")

      # get all the group normalized counts
      group_mean_expr <- rowData(norm_cts(myBbcEdgerobj)) %>%
        as.data.frame(stringsAsFactors = FALSE) %>%
        dplyr::select(dplyr::ends_with("norm_log_cpm")) %>%
        tibble::rownames_to_column("ensembl_id")

      # merge everything
      de_res <- top_tags %>%
        dplyr::left_join(gene_symbols, by = "ensembl_id") %>%
        dplyr::left_join(group_mean_expr, by = "ensembl_id")

      # sort by FDR and arrange the columns
      de_res <- de_res %>%
        dplyr::arrange(.data$FDR) %>%
        dplyr::select("gene_symbol", "ensembl_id",
                      dplyr::one_of(top_tags_colnames),
                      dplyr::ends_with("norm_log_cpm"))

      de_res
    })

    names(res_list) <- my_edger_res_names
  }


  res_list
}
