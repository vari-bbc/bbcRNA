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
#' @param BMdataset character. The name of the BM dataset you want to use. Try to
#'   use 'bbcRNA::searchBM("species_name")' to find the appropriate dataset if
#'   needed.
#' @return A BbcSE object.
#' @export
#' @importFrom AnnotationDbi keys mapIds
#' @importFrom biomaRt getBM useMart useDataset listDatasets
ens2entrez <- function(x, orgdb, BMdataset="") {
  if(!is(x, "BbcSE")) stop("x is not a BbcSE object")

  if(BMdataset==""){
    message("Using OrgDb.")
    if(!"ENSEMBL" %in% keytypes(orgdb)){
      stop("This org.db does not have Ensembl IDs. Use biomart (set BMdataset option).")
    }
    # get Ensembl IDs present in the OrgDb
    ens_genes <-  AnnotationDbi::keys(orgdb, keytype="ENSEMBL")

    # keep only genes present in BbcSE object
    ens_genes <- ens_genes[ens_genes %in% rownames(x)]

    # get entrez ids
    entrez_ids <- AnnotationDbi::mapIds(orgdb,
                                        keys = ens_genes,
                                        column = "ENTREZID",
                                        keytype = "ENSEMBL",
                                        multiVals = "first")

    stopifnot(identical(length(ens_genes), length(entrez_ids)))

    #names(entrez_ids) <- ens_genes

  } else {
    # select mart to use
    ensembl <- biomaRt::useMart("ensembl")

    # make sure dataset exists and print out version
    ensembl_datasets <- biomaRt::listDatasets(mart=ensembl)
    ensembl_dataset_to_use <- ensembl_datasets[ensembl_datasets$dataset==BMdataset, ]
    if(identical(nrow(ensembl_dataset_to_use), 0)){
      stop("BM dataset does not exist.")
    } else{
      stopifnot(nrow(ensembl_dataset_to_use) == 1) # make sure only one matching dataset
    }
    message(paste0("Using Biomart. Dataset version: ", ensembl_dataset_to_use$version))

    # get the dataset from BM
    ensembl <- biomaRt::useDataset(BMdataset, mart=ensembl)

    #ens_2_symbol <- biomaRt::getBM(attributes = c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
    ens_2_entrez <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = ensembl)

    # keep only genes present in BbcSE object
    ens_2_entrez <- ens_2_entrez[ens_2_entrez$ensembl_gene_id %in% rownames(x), ]

    # take the first match if there are multiple Entrez matches for a particular Ensembl ID
    ens_2_entrez <- ens_2_entrez[!duplicated(ens_2_entrez$ensembl_gene_id), ]

    entrez_ids <- ens_2_entrez$entrezgene_id
    names(entrez_ids) <- ens_2_entrez$ensembl_gene_id
  }

  message(paste0(sum(!rownames(x) %in% names(entrez_ids)),
                 " Ensembl IDs not in the database."))

  # remove Entrez genes that match more than one Ensembl gene
  dup_genes <- unique(entrez_ids[duplicated(entrez_ids)])

  message(paste0(sum(entrez_ids %in% dup_genes),
                 " Ensembl IDs matching the same Entrez ID as another Ensembl ID."))

  entrez_ids <- entrez_ids[!entrez_ids %in% dup_genes]

  # make sure entrez ids unique
  stopifnot(length(entrez_ids) ==  length(unique(entrez_ids)))

  message(paste0(length(entrez_ids),
                 " genes with valid Entrez ID out of ",
                 nrow(x),
                 " total genes"))

  # genes in the BbcSE object absent from OrgDb/Biomart or removed due to more
  # than one Ensembl match
  missing_genes <- rownames(x)[!rownames(x) %in% names(entrez_ids)]
  names(missing_genes) <- missing_genes # store the ensembl ID as names
  missing_genes[1:length(missing_genes)] <- NA # convert the ensembl IDs to NAs

  #combine the genes present in OrgDb with those absent
  entrez_ids <- c(entrez_ids, missing_genes)

  # order based on the genes in the BbcSE object
  rowData(x) <- cbind(rowData(x), entrez_ids = entrez_ids[rownames(x)])

  validObject(x)

  return(x)

}

###-----------------------------------------------------------------------------

#' Search biomart datasets with keyword
#'
#' Search biomart datasets with keyword
#'
#' @param regex regex to search for. Perl style.
#' @export
#' @importFrom biomaRt listDatasets
searchBM <- function(keyword){
  ensembl <- useMart("ensembl")

  bm_datasets <- biomaRt::listDatasets(ensembl)

  cols_w_keyword <- grep(keyword, bm_datasets, perl = TRUE, ignore.case = TRUE)
  rows_w_keyword_list <- lapply(cols_w_keyword, function(x){
    grep(keyword, bm_datasets[[x]], perl = TRUE, ignore.case = TRUE)
  })
  rows_w_keyword <- unique(do.call(c, rows_w_keyword_list))
  if(length(rows_w_keyword) == 0){
    stop("Keyword not found.")
  }
  return(bm_datasets[rows_w_keyword, ])
}

