#' Getter for aln_rates slot
#' @param x BbcSE object
#' @param withDimnames logical indicating whether output should have dim names
#' @export
aln_rates <- function(x, withDimnames=TRUE) {
  out <- x@aln_rates
  if (withDimnames) {
    rownames(out) <- colnames(x)
  }
  out
}

###-----------------------------------------------------------------------------

#' Setter for aln_rates slot
#'
#' @param x A BbcSE object.
#' @param value A matrix of mapping metrics. Rownames must correspond to samples
#' in counts matrix and column names indicating metric type must be set.
#' @export
`aln_rates<-` <- function(x, value) {

  # check that sample names match
  stopifnot(identical(nrow(value), ncol(x)))
  stopifnot(all(rownames(value) %in% colnames(x)))

  valid_aln_rates_colnames <- c("input_reads",
                                "uniq_aln_reads",
                                "mult_aln_reads")

  if (!all(valid_aln_rates_colnames %in% colnames(value)) &&
      identical(length(valid_aln_rates_colnames), length(colnames(value)))){
    stop(paste0("colnames(value) must be ",
                paste(valid_aln_rates_colnames, collapse = ", ")))
  }

  map_rate <-
    round(
      (value[, "uniq_aln_reads"] + value[, "mult_aln_reads"]) /
        value[, "input_reads"],
      digits = 3
    ) * 100

  uniq_map_rate <- round(value[, "uniq_aln_reads"] /
                           value[, "input_reads"],
                         digits = 3) * 100

  value <- cbind(value, map_rate, uniq_map_rate)

  # set the aln_rates slot. Order the aln_rates according to colnames(x)
  x@aln_rates <- value[colnames(x), ]
  validObject(x)
  x
}

###-----------------------------------------------------------------------------

#' Plot alignment metrics
#'
#' Plot alignment metrics as a bar chart.
#'
#' @param x A BbcSE object.
#' @param type aln_rates column to plot.
#' @param fill_by Variable to fill by. Can be column from colData() or
#'   aln_data().
#' @param facet_by Variable to facet by. Can be column from colData() or
#'   aln_data().
#' @return A ggplot object
#' @export
#' @import ggplot2
#' @importFrom tibble as_tibble
#' @importFrom cowplot theme_cowplot
plot_aln_rates <- function(x, type = "uniq_map_rate", fill_by = "sample",
                           facet_by = "") {
  if(!is(x, "BbcSE")) stop("x not BbcSE object")
  validObject(x)

  aln_data <- tibble::as_tibble(aln_rates(x),
                                rownames = "sample")

  sample_meta <- tibble::as_tibble(colData(x),
                                   rownames = "sample")

  if(!type %in% colnames(aln_data)){
    stop(paste0(type, " is not a column in aln_data"))
  }

  if(!fill_by %in% c(colnames(aln_data), colnames(sample_meta))){
    stop(paste0(type, " is not a column in aln_data or colData"))
  }

  if(!facet_by %in% c(colnames(aln_data), colnames(sample_meta))){
    stop(paste0(type, " is not a column in aln_data or colData"))
  }


  # join with sample meta data
  aln_data <- dplyr::left_join(aln_data, sample_meta, by = "sample")

  # divide read counts by 1M
  aln_data <-  aln_data %>%
    dplyr::mutate(input_reads = round(.data$input_reads / 10^6, digits=2),
                  uniq_aln_reads = round(.data$uniq_aln_reads / 10^6, digits=2),
                  mult_aln_reads = round(.data$mult_aln_reads / 10^6, digits=2))

  sensical_labels <- c(input_reads = "Total input reads (M)",
                       uniq_aln_reads = "Uniquely aligned reads (M)",
                       mult_aln_reads = "Multi aligned reads (M)",
                       map_rate = "Overall mapping rate",
                       uniq_map_rate = "Unique mapping rate")

  myplot <- ggplot(aln_data, aes_string(x="sample", y=type)) +
    geom_col(aes_string(fill=fill_by), color="black") +
    geom_text(aes_string(label=type), angle=45, vjust=-0.5, hjust=0) +
    ylab(sensical_labels[type]) +
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))

  if(type %in% c("map_rate", "uniq_map_rate")) {
    myplot <- myplot + ylim(0,105)
  }

  if(!identical(facet_by, "")){
    myplot <- myplot + facet_wrap(facet_by, nrow = 1,
                                  scales="free_x")
  }

  myplot
}

###-----------------------------------------------------------------------------


