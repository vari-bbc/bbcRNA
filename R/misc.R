#' PCA, screeplot and PERMANOVA test
#'
#' Runs PCA on normalized counts. Makes PCA and scree plots. Tests for
#' significant clustering based on specific variables using PERMANOVA and
#' Euclidean distances.
#'
#' @param x A BbcSE object.
#' @param norm_cts_type "edger" or "deseq2"
#' @param assay_name Name of assay from the 'norm_cts' slot. Option was
#'   implemented with the idea of being able to plot batch-corrected data.
#' @param color_by colData column to color by for PCA plot
#' @param shape_by colData column to shape by for PCA plot
#' @param adonis logical for whether vegan::adonis should be run. If TRUE,
#'   adonis() will be run with Euclidean distance calculated from the same
#'   normalized counts as used for PCA.
#' @param adonis_by colData column to test using PERMANOVA. May be a vector of
#'   values; in this case, each variable is tested sequentially using adonis().
#' @return A list containing ggplot objects for 1. PCA 2. scree plot 3.
#'   prcomp()$x merged with meta data. Handy if visualizations of >PC2 are
#'   needed 4. a list of output from vegan::adonis()
#' @seealso \code{\link[vegan]{adonis}}
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom cowplot theme_cowplot
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @importFrom vegan adonis vegdist
#' @export
plot_PCA <- function(x, norm_cts_type = "edger", assay_name="norm_log_cpm",
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
    norm_counts <- assay(norm_cts(edger(x)), assay_name)
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
    scale_x_continuous(breaks = seq(1, 100, 2)) +
    xlab("Principal Component") +
    ylab("Prop. of Variance") +
    ggtitle("PCA Plot of Expression Profiling") +
    theme_cowplot() +
    theme(axis.title.y = element_text(vjust=1),
          plot.margin = unit(c(0,0,0,6), "mm"))

  all_pca_output <- list(pca_plot, var_plot, pr_comps)

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
    all_pca_output <- list(pca_plot, var_plot, pr_comps, adonis_out)
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
#' @param assay_name Name of assay from the 'norm_cts' slot. Not compatible with
#'   grouped=TRUE. Option was implemented with the idea of being able to plot
#'   batch-corrected data.
#' @param gene_labels a column name from rowData. Must contain unique values
#' @param coldata_annot character vector of colData colnames for annotating
#' @param coldata_split character value of colData colname for splitting
#' @param rowdata_annot character vector of rowData colnames for annotating
#' @param rowdata_split character value of rowData colname for splitting
#' @param clust_rows logical for whether rows should be clustered
#' @param clust_cols logical for whether columns should be clustered
#' @param clustering_distance_rows see 'clustering_distance_rows' in \code{\link[ComplexHeatmap]{ComplexHeatmap::Heatmap}}
#' @param clustering_distance_columns see 'clustering_distance_columns' in \code{\link[ComplexHeatmap]{ComplexHeatmap::Heatmap}}
#' @param grouped logical indicating whether the mean expression per group
#'   should be shown.
#' @param zscores logical indicating whether the expression matrix should be
#'   converted to Z-scores
#' @return A Heatmap-class object.
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#'   \href{https://jokergoo.github.io/ComplexHeatmap-reference/book/}{ComplexHeatmap
#'    book}
#' @importFrom dplyr select ends_with
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation
#' @export
plot_heatmap <- function(x,
                         genes = rownames(x),
                         de_method = "edger",
                         assay_name = "norm_log_cpm",
                         gene_labels = NULL,
                         coldata_annot = NULL,
                         coldata_split = NULL,
                         rowdata_annot = NULL,
                         rowdata_split = NULL,
                         clust_rows = TRUE,
                         clust_cols = TRUE,
                         clustering_distance_rows = "euclidean",
                         clustering_distance_columns = "euclidean",
                         zscores = FALSE,
                         grouped = FALSE) {
  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  if(identical("edger", de_method)){
    curr_bbcedger <- edger(x)
    if(identical(BbcEdgeR(), curr_bbcedger)) stop("Empty 'edger' slot")

    # check if there are normalized counts for all requested genes
    genes_avail <- genes %in% rownames(dgelist(curr_bbcedger))
    if(!all(genes_avail)){
      if(sum(genes_avail) == 0) stop("No genes with normalized counts available.")
      num_genes_no_norm_cts <- sum(!genes_avail)
      message(paste0("Not plotting ", sum(!genes_avail),
                     " genes because normalized counts not available."))
      message(paste0("Plotting ", sum(genes_avail),
                     " genes with normalized counts."))
      genes <- genes[genes_avail]
    }

    # get the expression matrix
    if(isTRUE(grouped)){
      if(!is.null(coldata_annot) | !is.null(coldata_split)) {
        stop("Column data options not supported when grouped=TRUE")
      }
      if(assay_name != "norm_log_cpm") {
        stop("Only default 'assay_name' supported for grouped=TRUE")
      }
      expr_mat <- rowData(norm_cts(curr_bbcedger))[genes, , drop=FALSE] %>%
        as.data.frame(stringsAsFactors=FALSE) %>%
        dplyr::select(dplyr::ends_with(".norm_log_cpm")) %>%
        as.matrix()

    } else {
      expr_mat <- assay(norm_cts(curr_bbcedger), assay_name)[genes, , drop=FALSE]
    }

    # convert to zscores (gene-wise) if needed
    if(isTRUE(zscores)){
      expr_mat <- t(scale(t(expr_mat), center = TRUE, scale = TRUE))
      expr_scale_title <- "Z-scores"
    } else{
      expr_scale_title <- "Norm. log CPM"
    }

    # get gene labels, either from rowData or from rownames of expr_mat
    if(!is.null(gene_labels)){
      if(length(gene_labels) > 1) stop("gene_labels can only be one column")
      if(gene_labels == "rownames"){
        gene_labels <- rownames(expr_mat)
        names(gene_labels) <- gene_labels
      } else {
        if(!gene_labels %in% colnames(rowData(x)))
          stop(paste0(gene_labels, " does not exist in rowData"))
        gene_labels <- rowData(x)[genes, gene_labels]
        names(gene_labels) <- genes
      }
      showGeneNames <- TRUE
    } else{
      showGeneNames <- FALSE
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
    expr_mat <- expr_mat[genes, , drop=FALSE]

    if(!is.null(gene_labels)){
      gene_labels <- gene_labels[rownames(expr_mat)]
    }

    if(!is.null(coldata_split))
      coldata_split <- coldata_split[colnames(expr_mat)]

    if(!is.null(rowdata_split))
      rowdata_split <- rowdata_split[rownames(expr_mat)]

    # function that returns a list of named vectors which contain either colors
    # or a function for colors depending on whether the column is discrete or
    # continuous (numeric) respectively.
    prep_colors_for_complexheatmap <- function(df){
      df_isNumeric <- sapply(df, is.numeric)
      df_isFactor <- sapply(df, is.factor)
      rcolorbrwer_discrete_pals <- c("Dark2", "Paired", "Set1", "Accent", "Set2")

      color_list <- lapply(seq(1, length(df_isNumeric)), function(x){
        if(isFALSE(df_isNumeric[x])){
          if(isFALSE(df_isFactor[x])){
            uniq_values <- unique(df[, x])
          } else{
            uniq_values <- levels(df[, x])
          }

          if(length(uniq_values) > 8) stop("Annotation variable cannot have >8 unique values.")

          colors <- RColorBrewer::brewer.pal(8,
                                             rcolorbrwer_discrete_pals[x])[1:length(uniq_values)]
          names(colors) <- uniq_values
        } else{
          colors <-  circlize::colorRamp2(seq(min(df[, x]),
                                              max(df[, x]), length = 3),
                                          c("blue", "#EEEEEE", "red"))
        }
        return(colors)
      })

      names(color_list) <- colnames(df)

      color_list
    }
    ### End function

    if(!is.null(coldata_annot)){
      # complexheatmap assigns random colors to discrete annotations
      # here we define the palette manually
      coldata_annot <- ComplexHeatmap::HeatmapAnnotation(
        df = coldata_annot[colnames(expr_mat), , drop=FALSE],
        col = prep_colors_for_complexheatmap(coldata_annot))
    }

    if(!is.null(rowdata_annot)){
      # complexheatmap assigns random colors to discrete annotations
      # here we define the palette manually
      rowdata_annot <- ComplexHeatmap::rowAnnotation(
        df = rowdata_annot[rownames(expr_mat), , drop=FALSE],
        col = prep_colors_for_complexheatmap(rowdata_annot))
    }

    # make main heatmap
    expr_ht <- ComplexHeatmap::Heatmap(expr_mat,
                                       name = expr_scale_title,
                                       cluster_row_slices = FALSE,
                                       cluster_column_slices = FALSE,
                                       cluster_rows = clust_rows,
                                       cluster_columns = clust_cols,
                                       clustering_distance_rows = clustering_distance_rows,
                                       clustering_distance_columns = clustering_distance_columns,
                                       row_split = rowdata_split,
                                       column_split = coldata_split,
                                       top_annotation = coldata_annot,
                                       right_annotation = rowdata_annot,
                                       row_title_rot = 0,
                                       row_labels = gene_labels,
                                       show_row_names = showGeneNames)


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
    ggplot2::geom_histogram(color="black", fill="gray55",
                            breaks = seq(0, 1, 0.05)) +
    ggplot2::facet_wrap(c("contrast_name")) +
    cowplot::theme_cowplot()

  return(pval_plot)
}


###-----------------------------------------------------------------------------
#' Batch correct normalized counts
#'
#' Batch correct normalized counts for visualizations or other downstream non-DE
#' analysis applications.
#'
#' @param x A BbcSE object
#' @param de_method "edger" or "deseq2"
#' @param new_assay_name Name of the batch-corrected counts stored in new assay
#'   in 'norm_cts' slot.
#' @param correction_method "removeBatchEffect" or "combat"
#' @param ... passed to batch correction function.
#' @return A BbcSE object
#' @importFrom limma removeBatchEffect
#' @examples
#' \dontrun{
#' # Default is limma::removeBatchEffect
#' bbc_obj_batch <- batch_correct_norm_cts(bbc_obj, batch=colData(bbc_obj)$Rep,
#' design=model.matrix(~Condition, data=colData(bbc_obj)))
#'
#' plot_PCA(bbc_obj_batch, assay_name = "batch_corr", adonis=FALSE,
#' color_by="Time", shape_by="Rep")
#'
#' # Combat is also supported
#' bbc_obj_batch <- batch_correct_norm_cts(bbc_obj, new_assay_name="combat",
#' batch=colData(bbc_obj)$Rep, correction_method = "combat", mod =
#' model.matrix(~Condition, data=colData(bbc_obj)))
#'
#' plot_PCA(bbc_obj_batch, assay_name = "combat", adonis=FALSE, color_by="Time",
#' shape_by="Rep")
#' }
#' @seealso \code{\link[limma]{removeBatchEffect} \link[sva]{ComBat}}
#' @export
batch_correct_norm_cts <- function(x,
                              de_method = "edger",
                              correction_method = "removeBatchEffect",
                              new_assay_name = "batch_corr",
                              ...) {

  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  if (identical("combat", correction_method) &
      !requireNamespace("sva", quietly = TRUE)) {
    stop("Package \"sva\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if(identical("edger", de_method)){
    # should never be triggered but this is a fail-safe that allows colData()
    # columns to be safely used as batch factors.
    if(!identical(rownames(colData(x)), colnames(x@edger@norm_cts))){
      stop("Samples in BbcSE meta data do not match those in 'norm_cts' slot.")
    }

    emat <- assay(x@edger@norm_cts, "norm_log_cpm")

    if(identical("removeBatchEffect", correction_method)){
      if(length(list(...)) == 0) stop("See ?limma::removeBatchEffect for arguments to perform batch correction.")
      assay(x@edger@norm_cts, new_assay_name) <-
        limma::removeBatchEffect(emat, ...)
    } else if (identical("combat", correction_method)){

      if(length(list(...)) == 0) stop("See ?sva::ComBat for arguments to perform batch correction.")
      assay(x@edger@norm_cts, new_assay_name) <- sva::ComBat(dat=emat, ...)
    }
    else{
      stop("Unsupported 'correction_method'.")
    }
  } else{
    stop("Invalid 'de_method'.")
  }

  return(x)
}

###-----------------------------------------------------------------------------
#' Volcano plot
#'
#' Plot volcano plots for contrasts in edger or deseq2 slot of BbcSE object
#'
#' @param x A BbcSE object
#' @param de_method "edger" or "deseq2"
#' @param gene_name "gene_symbol" or "ensembl_id". Corresponds to columns in
#'   'bbcRNA::get_de_table()' output
#' @param pval_name "PValue" or "FDR". Corresponds to columns in
#'  'bbcRNA::get_de_table()' output
#' @param contrast_names character value or vector for the contrast(s) of
#'   interest. See name(de_results(edger(x))). If not specified, then all
#'   contrasts will be processed.
#' @param pvalCutoff "Cut-off for statistical significance. A horizontal line
#'   will be drawn at -log10(pCutoff)." See EnhancedVolcano documentation.
#' @param y_title Y axis label
#' @param ... Passed to EnhancedVolcano::EnhancedVolcano
#' @return A list of ggplot objects
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @seealso \code{\link[EnhancedVolcano]{EnhancedVolcano}}
#' @export
plot_volcano <- function(x,
                         de_method = "edger",
                         gene_name = "gene_symbol",
                         pval_name = "FDR",
                         pvalCutoff = 0.05,
                         y_title = paste0("-Log10 (", pval_name, ")"),
                         contrast_names,
                         ...) {

  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  de_table_list <- get_de_table(x, de_pkg = de_method)

  # make the gene IDs the row names
  de_table_list <- lapply(de_table_list, function(table){
    rownames(table) <- table[, gene_name]
    return(table)
  })

  if(identical("edger", de_method)){

    if (!missing(contrast_names)){
      edger_contrast_names <- contrast_names
      if(!all(edger_contrast_names %in% names(de_table_list))){
        stop("Check that 'contrast_names' all exist in the BbcSE object.")
      }
    } else{
      edger_contrast_names <- names(de_table_list)
    }

    volcano_plots <- lapply(edger_contrast_names, function(table){
      EnhancedVolcano::EnhancedVolcano(toptable=de_table_list[[table]],
                                       lab=rownames(de_table_list[[table]]),
                                       x="logFC",
                                       y=pval_name,
                                       title = table,
                                       subtitle = "",
                                       ylab=y_title,
                                       ylim = c(0, max(max(-log10(de_table_list[[table]][,pval_name]),
                                                           na.rm=TRUE),
                                                           -log10(pvalCutoff))),
                                       subtitleLabSize = 0,
                                       caption = paste0('Total = ',
                                                        nrow(de_table_list[[table]]),
                                                        ' genes'),
                                       pCutoff=pvalCutoff, ...)
    })
    names(volcano_plots) <- edger_contrast_names

  }

  return(volcano_plots)
}
