#' Convert Ensembl IDs to gene symbols.
#'
#' Convert Ensembl IDs to gene symbols.
#'
#' Gene symbols with more than one Ensembl ID match are determined only from
#' genes present in the BbcSE object, so be sure it is unfiltered.
#'
#' \itemize{
#'   \item Symbols that match multiple Ensembl IDs will be resolved by concatenating
#' the Ensembl ID and the gene symbol. Uses scater::uniquifyFeatureNames.
#'   \item Genes with multiple possible symbols will be labelled as NA.
#'   \item Genes absent from OrgDb will be labelled as Ensembl ID.
#' }
#'
#' @param x A BbcSE object.
#' @param orgdb A OrgDb object. Download the annotations for your species from
#'   'http://bioconductor.org/packages/release/BiocViews.html#___OrgDb'
#' @return A BbcSE object.
#' @export
#' @importFrom AnnotationDbi keys mapIds
ens2sym <- function(x, orgdb) {
  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  # get Ensembl IDs present in the OrgDb
  ens_genes <-  AnnotationDbi::keys(orgdb, keytype="ENSEMBL")

  # keep only genes present in BbcSE object
  ens_genes <- ens_genes[ens_genes %in% rownames(x)]

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

  uniq_syms <- uniquifyFeatureNames(names(gene_syms), gene_syms)

  names(uniq_syms) <- names(gene_syms)

  # genes in the BbcSE object absent from OrgDb
  missing_syms <- rownames(x)[!rownames(x) %in% names(uniq_syms)]
  names(missing_syms) <- missing_syms

  #combine the genes present in OrgDb with those absent
  uniq_syms <- c(uniq_syms, missing_syms)

  # keep only the genes present in the BbcSE object and order based on
  # the genes in the BbcSE object
  rowData(x) <- cbind(rowData(x), uniq_syms = uniq_syms[rownames(x)])

  validObject(x)

  return(x)

}

###-----------------------------------------------------------------------------

#' Convert Ensembl IDs to Entrez IDs for gene enrichment analyses.
#'
#' Convert Ensembl IDs to Entrez IDs for gene enrichment analyses.
#'
#' Entrez IDs with more than one Ensembl ID match are determined only from genes
#' present in the BbcSE object, so be sure it is unfiltered.
#'
#' \itemize{
#' \item Entrez IDs that match multiple Ensembl IDs are NA
#' \item The first Entrez ID will be returned for genes with multiple possible
#' Entrez IDs.
#' \item Genes with no Entrez ID or absent from the database will be
#' NA. }
#'
#' @param x A BbcSE object.
#' @param orgdb A OrgDb object. Download the annotations for your species from
#'   'http://bioconductor.org/packages/release/BiocViews.html#___OrgDb'
#' @return A BbcSE object.
#' @export
#' @importFrom AnnotationDbi keys mapIds
ens2entrez <- function(x, orgdb) {
  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  # get Ensembl IDs present in the OrgDb
  ens_genes <-  AnnotationDbi::keys(orgdb, keytype="ENSEMBL")

  # keep only genes present in BbcSE object
  ens_genes <- ens_genes[ens_genes %in% rownames(x)]

  entrez_ids <- AnnotationDbi::mapIds(orgdb,
                                     keys = ens_genes,
                                     column = "ENTREZID",
                                     keytype = "ENSEMBL",
                                     multiVals = "first")

  stopifnot(identical(length(ens_genes), length(entrez_ids)))

  names(entrez_ids) <- ens_genes

  # remove Entrez genes that match more than one Ensembl gene
  dup_genes <- unique(entrez_ids[duplicated(entrez_ids)])
  entrez_ids <- entrez_ids[!entrez_ids %in% dup_genes]

  message(paste0(length(entrez_ids),
                 " genes with valid Entrez ID out of ",
                 nrow(x),
                 " total genes"))

  # genes in the BbcSE object absent from OrgDb or removed due to more than one
  # Ensembl match
  missing_genes <- rownames(x)[!rownames(x) %in% names(entrez_ids)]
  names(missing_genes) <- missing_genes # store the ensembl ID as names
  missing_genes <- NA # convert the ensembl IDs to NAs

  #combine the genes present in OrgDb with those absent
  entrez_ids <- c(entrez_ids, missing_genes)

  # keep only the genes present in the BbcSE object and order based on
  # the genes in the BbcSE object
  rowData(x) <- cbind(rowData(x), entrez_ids = entrez_ids[rownames(x)])

  validObject(x)

  return(x)

}

