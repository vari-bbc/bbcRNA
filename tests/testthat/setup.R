library(bbcRNA)
library(SummarizedExperiment)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(stringr)
#library(org.Ss.eg.db) # for testing biomart functionality for ens2entrz and ens2sym

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
aln_metrics <- read_star_aln_metrics(dir = ext_data_dir,
                                     rgx = "^[^_]+_[^_]+")

# get column meta data
col_meta <- read_col_meta(paste0(ext_data_dir, "/meta.txt"))
col_meta$genotype <- factor(col_meta$genotype, levels=c("WT", "mut"))

# make BbcSE object
bbc_obj <- BbcSE(counts = counts_mat, aln_metrics = aln_metrics,
                 colData = col_meta[colnames(counts_mat), ])

bbc_obj_no_aln_metrics <- BbcSE(counts = counts_mat,
                                colData = col_meta[colnames(counts_mat), ])

# filter lowly expressed genes, calculate norm factors, make DGEList
bbc_obj_dgelist <- makeDGEList(bbc_obj, group="genotype")

# Identify DE genes with glmQLFTest
bbc_obj_glmQLFTest <- findDEGs(bbc_obj_dgelist,
                               de_pkg = "edger",
                               test = "glmQLFTest",
                               design = "~0+genotype",
                               contrasts = list(c("genotype", "mut", "WT")))
bbc_obj_glmQLFTest <- ens2sym(bbc_obj_glmQLFTest, org.Mm.eg.db)
bbc_obj_glmQLFTest <- ens2entrez(bbc_obj_glmQLFTest, org.Mm.eg.db)


# Identify DE genes with glmTreat
bbc_obj_glmTreat <- findDEGs(bbc_obj_dgelist,
                             de_pkg = "edger",
                             test = "glmTreat",
                             design = "~0+genotype",
                             contrasts = list(c("genotype", "mut", "WT")))

bbc_obj_glmTreat <- ens2sym(bbc_obj_glmTreat, org.Mm.eg.db)

# Identify DE genes using coefs
bbc_obj_glmQLFTest_coef <- findDEGs(x=bbc_obj_dgelist, design="~genotype",
                                    coefs=list(2))

bbc_obj_glmTreat_coef <- findDEGs(x=bbc_obj_dgelist, design="~genotype",
                                  coefs=list(2), test="glmTreat", lfc = log2(2))


