#' PCA, screeplot and PERMANOVA test
#'
#' Runs PCA on normalized counts. Makes PCA and scree plots. Tests for
#' significant clustering based on specific variables using PERMANOVA and
#' Euclidean distances.
#'
#' @param x A BbcSE object.
#' @param norm_cts_type "edger" or "deseq2"
#' @param color_by colData column to color by for PCA plot
#' @param shape_by colData column to shape by for PCA plot
#' @param adonis logical for whether vegan::adonis should be run. If TRUE,
#'   adonis() will be run with Euclidean distance calculated from the same
#'   normalized counts as used for PCA.
#' @param adonis_by colData column to test using PERMANOVA. May be a vector of
#'   values; in this case, each variable is tested sequentially using adonis().
#' @return A list containing ggplot objects for 1. PCA 2. scree plot 3. a list
#'   of output from vegan::adonis()
#' @seealso \code{\link[vegan]{adonis}}
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom cowplot theme_cowplot
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom vegan adonis vegdist
#' @export
plot_PCA <- function(x, norm_cts_type = "edger",
                     color_by = colnames(colData(x))[[1]],
                     shape_by = colnames(colData(x))[[1]],
                     adonis = TRUE,
                     adonis_by = colnames(colData(x))[[1]]
                     ){
  if (!is(x, "BbcSE")) stop("x is not a BbcSE object")

  if (length(color_by) > 1) {
    stop("'color_by' must have length of 1")
  }
  if (length(shape_by) > 1) {
    stop("'shape_by' must have length of 1")
  }
  if (!color_by %in% colnames(colData(x))) {
    stop("'color_by' parameter not found in column meta data")
  }
  if (!shape_by %in% colnames(colData(x))) {
    stop("'shape_by' parameter not found in column meta data")
  }
  if (!all(adonis_by %in% colnames(colData(x)))) {
    stop("One or more 'adonis_by' parameter(s) not found in column meta data")
  }

  if (norm_cts_type == "edger"){
    norm_counts <- assay(norm_cts(edger(x)))
    pca <- prcomp(t(norm_counts))
  } else if (norm_cts_type == "deseq2"){
    # not implemented yet
  }

  pr_comps <- data.frame(pca$x)
  pr_comps$sample <- rownames(pr_comps)

  column_meta <- as.data.frame(colData(x), stringsAsFactor=FALSE) %>%
    tibble::rownames_to_column("sample")

  pr_comps <- dplyr::left_join(pr_comps, column_meta, by="sample")

  # Plot percent variation explained
  prop_var <- data.frame(t(summary(pca)$importance))
  names(prop_var) = c('sd', 'prop', 'cum')
  prop_var$num = 1:nrow(prop_var)

  pca_plot <- ggplot(pr_comps, aes_string(x="PC1", y="PC2")) + #
    geom_point(size=3, aes_string(color=color_by, pch=shape_by)) +
    xlab(paste0("PC1 (", prop_var[prop_var$num == 1, "prop"]*100, "%)")) +
    ylab(paste0("PC2 (", prop_var[prop_var$num == 2, "prop"]*100, "%)")) +
    theme_cowplot() +
    scale_color_brewer(palette = 'Paired')

  var_plot <- ggplot(prop_var %>% dplyr::filter(.data$num <= 12),
                     aes_string(x="num", y="prop")) +
    geom_point(size=1.5) +
    geom_line() +
    xlab("Principal Component") +
    ylab("Prop. of Variance") +
    ggtitle("PCA Plot of Expression Profiling") +
    theme_cowplot() +
    theme(axis.title.y = element_text(vjust=1),
          plot.margin = unit(c(0,0,0,6), "mm"))

  all_pca_output <- list(pca_plot, var_plot)

  if (isTRUE(adonis)){
    # calculate Euclidean distance
    eucl_dist <- vegan::vegdist(t(norm_counts), method="euclidean")
    column_meta <- colData(x)

    # run adonis() for each 'adonis_by' variable. I couldn't figure out how to
    # specify the RHS of the 'formula' paramter for vegan::adonis, so as a
    # work-around, I rename the desired variable to 'group' for each lapply
    # iteration.
    adonis_out <- lapply(adonis_by, function(curr_col){
      column_meta_temp <- data.frame(group = column_meta[[curr_col]],
                                     stringsAsFactors = FALSE)

      vegan::adonis(eucl_dist ~ group,
                    data = column_meta_temp, method = "eu")
    })
    names(adonis_out) <-  adonis_by

    # get the P values
    adonis_pvals <- sapply(adonis_by, function(curr_col) {
      round(adonis_out[[curr_col]]$aov.tab$`Pr(>F)`[[1]], digits = 3)
    })
    names(adonis_pvals) <- adonis_by

    # add the P values to the PCA plot as a caption
    pvals_concat <- paste0("P-values: ",
                           paste0(names(adonis_pvals), "=", adonis_pvals,
                           collapse = ", "))
    pca_plot <- pca_plot + labs(caption = pvals_concat)

    # assemble output
    all_pca_output <- list(pca_plot, var_plot, adonis_out)
  }

  all_pca_output
}


###-----------------------------------------------------------------------------
#' Plot heatmap
#'
#' Plot heatmap with option for splitting rows and columns and adding
#' annotation.
#'
#' @param x A BbcSE object
#' @param genes gene IDs (rownames)
#' @param de_method "edger" or "deseq2"
#' @param gene_labels a column name from rowData. Must contain unique values
#' @param coldata_annot character vector of colData colnames for annotating
#' @param coldata_split character value of colData colname for splitting
#' @param rowdata_annot character vector of rowData colnames for annotating
#' @param rowdata_split character value of rowData colname for splitting
#' @param clust_rows logical for whether rows should be clustered
#' @param clust_cols logical for whether columns should be clustered
#' @param grouped logical indicating whether the mean expression per group
#'   should be shown.
#' @param zscores logical indicating whether the expression matrix should be
#'   converted to Z-scores
#' @return A Heatmap-class list.
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#'   \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/}{ComplexHeatmap
#'    book}
#' @importFrom dplyr select ends_with
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @export
plot_heatmap <- function(x,
                         genes = rownames(x),
                         de_method = "edger",
                         gene_labels = NULL,
                         coldata_annot = NULL,
                         coldata_split = NULL,
                         rowdata_annot = NULL,
                         rowdata_split = NULL,
                         clust_rows = TRUE,
                         clust_cols = TRUE,
                         zscores = FALSE,
                         grouped = FALSE) {

  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  if(identical("edger", de_method)){
    curr_bbcedger <- edger(x)
    if(identical(BbcEdgeR(), curr_bbcedger)) stop("Empty 'edger' slot")

    # get the expression matrix
    if(isTRUE(grouped)){
      if(!is.null(coldata_annot) | !is.null(coldata_split)) {
       stop("column data options not compatible with grouped=TRUE")
      }
      expr_mat <- rowData(norm_cts(curr_bbcedger))[genes, ] %>%
        as.data.frame(stringsAsFactors=FALSE) %>%
        dplyr::select(dplyr::ends_with(".norm_log_cpm")) %>%
        as.matrix()

    } else {
      expr_mat <- assay(norm_cts(curr_bbcedger), "norm_log_cpm")[genes, ]
    }

    # convert to zscores (gene-wise) if needed
    if(isTRUE(zscores)){
      expr_mat <- t(scale(t(expr_mat), center = TRUE, scale = TRUE))
    }

    # get gene labels, either from rowData or from rownames of expr_mat
    if(!is.null(gene_labels)){
      if(length(gene_labels) > 1) stop("gene_labels can only be one column")
      gene_labels <- rowData(x)[genes, gene_labels]
    } else{
      gene_labels <- rownames(expr_mat)
      names(gene_labels) <- gene_labels
    }

    ## make sure no duplicated gene labels
    stopifnot(length(gene_labels) == length(unique(gene_labels)))

    # coldata_split
    if(!is.null(coldata_split)){
      if(length(coldata_split) > 1) stop("coldata_split can only be one column")
      if(!coldata_split %in% colnames(colData(x))) {
        stop("Requested coldata_split not in colData.")
      }
      coldata_split <- colData(x)[colnames(expr_mat), coldata_split,
                                  drop = TRUE]
      names(coldata_split) <- colnames(expr_mat)
    }

    # col annots
    if(!is.null(coldata_annot)){
      coldata_valid <- coldata_annot %in% colnames(colData(x))
      if(!all(coldata_valid)) {
        stop(paste("Requested coldata_annot not in colData:",
                   coldata_annot[!coldata_valid],
                   sep = " "))
      }
      coldata_annot <- colData(x)[colnames(expr_mat), coldata_annot,
                                  drop = FALSE] %>%
        as.data.frame(stringsAsFactors = FALSE)

    }

    # rowdata_split
    if(!is.null(rowdata_split)){
      if(length(rowdata_split) > 1) stop("rowdata_split can only be one row")
      if(!rowdata_split %in% colnames(rowData(x))) {
        stop("Requested rowdata_split not in rowData.")
      }
      rowdata_split <- rowData(x)[genes, rowdata_split, drop = TRUE]
      names(rowdata_split) <- genes
    }

    # row annots
    if(!is.null(rowdata_annot)){
      rowdata_valid <- rowdata_annot %in% colnames(rowData(x))
      if(!all(rowdata_valid)) {
        stop(paste("Requested rowdata_annot not in rowData:",
                   rowdata_annot[!rowdata_valid],
                   sep = " "))
      }
      rowdata_annot <- rowData(x)[rownames(expr_mat), rowdata_annot,
                                  drop = FALSE] %>%
        as.data.frame(stringsAsFactors = FALSE)

    }

    # double-check things are in order even though they should be sorted already
    expr_mat <- expr_mat[genes, ]
    gene_labels <- gene_labels[rownames(expr_mat)]

    if(!is.null(coldata_split))
      coldata_split <- coldata_split[colnames(expr_mat)]

    if(!is.null(rowdata_split))
      rowdata_split <- rowdata_split[rownames(expr_mat)]

    if(!is.null(coldata_annot)){
      coldata_annot <- ComplexHeatmap::HeatmapAnnotation(
        df = coldata_annot[colnames(expr_mat),])
    }

    if(!is.null(rowdata_annot)){
      rowdata_annot <- ComplexHeatmap::rowAnnotation(
        df = rowdata_annot[rownames(expr_mat),])
    }

    # make main heatmap
    expr_ht <- ComplexHeatmap::Heatmap(expr_mat,
                                       cluster_row_slices = FALSE,
                                       cluster_column_slices = FALSE,
                                       cluster_rows = clust_rows,
                                       cluster_columns = clust_cols,
                                       column_split = coldata_split,
                                       row_split = rowdata_split,
                                       top_annotation = coldata_annot,
                                       right_annotation = rowdata_annot,
                                       row_title_rot = 0,
                                       row_labels = gene_labels)


  } else if (identical("deseq2", de_method))
  {
    stop("deseq2 option not implemented yet")
  }

  return(expr_ht)
}

###-----------------------------------------------------------------------------
#' Plot P-values
#'
#' Plot P-value distribution for contrasts in edger or deseq2 slot of BbcSE
#' object
#'
#' @param x A BbcSE object
#' @param de_method "edger" or "deseq2"
#' @param contrast_names character value or vector for the contrast(s) of
#'   interest. See name(de_results(edger(x))). If not specified, then all
#'   contrasts will be processed.
#' @return A ggplot object
#' @import ggplot2
#' @importFrom cowplot theme_cowplot
#' @export
plot_pval_distrib <- function(x,
                         de_method = "edger",
                         contrast_names) {

  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  if(identical("edger", de_method)){
    # first element is for the model fit. Remove it.
    gene_and_pval <- de_results(edger(x))[-1]

    if (!missing(contrast_names)){
      edger_contrast_names <- contrast_names
      if(!all(edger_contrast_names %in% names(gene_and_pval))){
       stop("Check that 'contrast_names' all exist in the BbcSE object.")
      }
    } else{
      edger_contrast_names <- names(gene_and_pval)
    }

    # extract the P-values and add column for the contrast name
    pvals <- lapply(edger_contrast_names, function(contrast_name){
      out_df <- gene_and_pval[[contrast_name]]$table[, "PValue", drop = FALSE]
      out_df$contrast_name <- contrast_name
      out_df
    })
    pvals_combined <- do.call(rbind, pvals)

  } else if (identical("deseq2", de_method))
  {
    stop("deseq2 option not implemented yet")
    gene_and_pval <- ""
  }

  pvals_combined <- pvals_combined %>%
    dplyr::group_by(.data$contrast_name) %>%
    dplyr::mutate(gene_ct = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(contrast_name = paste0(.data$contrast_name,
                                         " (n=",
                                         format(.data$gene_ct, big.mark=",",
                                                scientific=FALSE),
                                         ")"))

  pval_plot <- ggplot(data = pvals_combined,
                      ggplot2::aes_string(x = "PValue")) +
    ggplot2::geom_histogram(color="black", fill="gray55", bins = 20) +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1)) +
    ggplot2::facet_wrap(c("contrast_name")) +
    cowplot::theme_cowplot()

  return(pval_plot)
}
