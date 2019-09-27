library(bbcRNA)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)

# directory containing external data
ext_data_dir <- system.file("extdata/tcell", package = "bbcRNA")

# get path to the count files
count_files <- grep(".*ReadsPerGene.out.tab(\\.[^\\.\\s]+)?$",
                    list.files(ext_data_dir, full.names = TRUE),
                    value = TRUE, perl = TRUE)

# get path to the alignment metrics files
aln_metrics_files <- grep(".*Log.final.out(\\.[^\\.\\s]+)?$",
                          list.files(ext_data_dir, full.names = TRUE),
                          value = TRUE, perl = TRUE)

# get counts matrix
counts_mat <- star_to_mat(dir = ext_data_dir,
                          rgx = "^[^_]+_[^_]+", column = 2)

# get the alignment metrics
aln_rates <- read_star_map_rates(dir = ext_data_dir,
                                 rgx = "^[^_]+_[^_]+")

# get column meta data
col_meta <- read_col_meta(paste0(ext_data_dir, "/meta.txt"))

# make BbcSE object
bbc_obj <- BbcSE(counts = counts_mat)



