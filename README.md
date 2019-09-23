
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bbcRNA

The goal of bbcRNA is to facilitate RNA-seq analysis, starting from
counts and ending with differentially expressed gene lists and
enerichment analyses.

## Installation

You can install bbcRNA from [GitHub](https://github.com/vari-bbc) with:

``` r
devtools::install_github("vari-bbc/bbcRNA")
```

## The BbcSE class

The BbcSE class extends SummarizedExperiment. Additional features
include:

  - An aln\_rates slot to hold a matrix containing alignment metrics for
    samples.
  - An edger slot for storing DGEList and edgeR result objects
  - A deseq2 slot for storing DESeqDataSet and DESeqResults

<!-- end list -->

``` r
# attach the package
library(bbcRNA)

# run the BbcSE constructor without any data to show the structure of the object

BbcSE()
#> class: BbcSE 
#> dim: 0 0 
#> metadata(0):
#> assays(1): counts
#> rownames: NULL
#> rowData names(0):
#> colnames: NULL
#> colData names(0):
#> aln_rates(0): 
#> edger( 0 ):   
#> deseq2( 0 ):
```
