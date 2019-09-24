#' Make PCA and scree plots
#'
#' Make PCA and scree plots.
#'
#' @param x A BbcSE object.
#' @param norm_cts_type "edger" or "deseq2"
#' @param color_by colData column to color by for PCA plot
#' @param shape_by colData column to shape by for PCA plot
#' @return A list containing ggplot objects for 1. PCA 2. scree plot.
#' @import ggplot2
#' @importFrom stats prcomp
#' @importFrom cowplot theme_cowplot
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr left_join
#' @export
plot_PCA <- function(x, norm_cts_type = "edger", color_by, shape_by){
  if (missing(color_by)) stop("Provide color_by parameter")
  if (missing(shape_by)) stop("Provide shape_by parameter")
  if (norm_cts_type == "edger"){
    norm_counts <- assay(norm_cts(edger(x)))
    pca <- prcomp(t(norm_counts))
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

  all_pca_plots <- list(pca_plot, var_plot)

  all_pca_plots
}
