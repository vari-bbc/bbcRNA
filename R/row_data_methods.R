#' Convert Ensembl IDs to gene symbols.
#'
#' Convert Ensembl IDs to gene symbols.
#'
#' Non 1:1 matches are determined from all the Ensembl IDs in the OrgDb, not
#' just those in the BbcSE object.
#' \itemize{
#'   \item Symbols that match multiple Ensembl IDs will be resolved by concatenating
#' the Ensembl ID and the gene symbol. Uses scater::uniquifyFeatureNames.
#'   \item Genes with multiple possible symbols will be labelled as NA.
#'   \item Genes absent from OrgDb will be labelled as Ensembl ID.
#' }
#'
#' @name ens2sym
#' @param x A BbcSE object.
#' @param orgdb A OrgDb object. Download the annotations for your species from
#'   'http://bioconductor.org/packages/release/BiocViews.html#___OrgDb'
#' @importFrom AnnotationDbi keys mapIds
#' @return A BbcSE object.
setMethod("ens2sym", c(x = "BbcSE"),
          function(x, orgdb) {
            ens_genes <-  AnnotationDbi::keys(orgdb, keytype="ENSEMBL")

            gene_syms <- AnnotationDbi::mapIds(orgdb,
                                               keys = ens_genes,
                                               column = "SYMBOL",
                                               keytype = "ENSEMBL",
                                               multiVals = "asNA")

            stopifnot(identical(length(ens_genes), length(gene_syms)))

            # function to uniquify gene names from Scater
            uniquifyFeatureNames <- function (ID, names)
            {
              if (length(ID) != length(names)) {
                stop("lengths of 'ID' and 'names' must be equal")
              }
              missing.name <- is.na(names)
              names[missing.name] <- ID[missing.name]
              dup.name <- names %in% names[duplicated(names)]
              names[dup.name] <- paste0(names[dup.name], "_", ID[dup.name])
              return(names)
            }

            uniq_syms <- uniquifyFeatureNames(ens_genes, gene_syms)

            names(uniq_syms) <- ens_genes

            # genes in the BbcSE object absent from OrgDb
            missing_syms <- rownames(x)[!rownames(x) %in% ens_genes]
            names(missing_syms) <- missing_syms

            #combine the genes present in OrgDb with those absent
            uniq_syms <- c(uniq_syms, missing_syms)

            # keep only the genes present in the BbcSE object and order based on
            # the genes in the BbcSE object
            rowData(x) <- cbind(rowData(x), uniq_syms = uniq_syms[rownames(x)])

            return(x)

          }
)
