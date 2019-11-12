
test_that("ens2sym works", {
  test_gene_syms <- function(bbc, ens2syms){

    bbc_rowdata <- rowData(bbc)

    # add genes that are absent from database
    missing_genes <- rownames(bbc_rowdata)[!rownames(bbc_rowdata) %in%
                                             names(ens2syms)]
    names(missing_genes) <- missing_genes

    gene_syms_df <- data.frame(symbols = c(ens2syms,
                                           rep(NA, length(missing_genes))),
                               stringsAsFactors = FALSE)

    rownames(gene_syms_df) <- c(names(ens2syms), names(missing_genes))

    # replace NAs with Ensembl ID
    gene_syms_df$symbols <- ifelse(is.na(gene_syms_df$symbols),
                                   rownames(gene_syms_df), gene_syms_df$symbols)

    # concatenate gene symbols with Ensembl IDs
    gene_syms_df$concat_sym <- paste0(gene_syms_df$symbols, "_",
                                      rownames(gene_syms_df))

    # match gene order
    stopifnot(all(rownames(bbc_rowdata) %in% rownames(gene_syms_df)))
    gene_syms_df <- gene_syms_df[rownames(bbc_rowdata), ]

    # test non-uniquified symbols correct
    bbc_uniquified <- stringr::str_detect(bbc_rowdata$uniq_syms,
                                          "_ENS[A-Z]*G[0-9]+$")

    expect_equivalent(bbc_rowdata[!bbc_uniquified, "uniq_syms"],
                      gene_syms_df[!bbc_uniquified, "symbols"])

    # test uniquified symbols correct
    expect_equivalent(bbc_rowdata[bbc_uniquified, "uniq_syms"],
                      gene_syms_df[bbc_uniquified, "concat_sym"])

    # test uniquified symbols are duplicated
    expect_true(all(table(gene_syms_df[bbc_uniquified, "symbols"]) > 1))

    invisible()
  }

  # OrgDb
  ens_genes <-  AnnotationDbi::keys(org.Mm.eg.db, keytype="ENSEMBL")

  gene_syms <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                     keys = ens_genes,
                                     column = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVals = "asNA")

  test_gene_syms(ens2sym(bbc_obj, org.Mm.eg.db), gene_syms)

  # Biomart
  ensembl <- biomaRt::useMart("ensembl")
  ensembl <- biomaRt::useDataset("mmusculus_gene_ensembl", mart=ensembl)
  ens_2_sym <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                                "external_gene_name"),
                                 mart = ensembl)

  # Make NA if there are multiple possible symbols
  ## Find Ensembl IDs with multiple symbols
  ens_w_multi_sym <- unique(ens_2_sym[duplicated(ens_2_sym$ensembl_gene_id),
                                         "ensembl_gene_id", drop=TRUE])
  ## Remove Ensembl IDs with multiple symbols
  ens_2_sym <- ens_2_sym[!ens_2_sym$ensembl_gene_id %in% ens_w_multi_sym, ]

  gene_syms_BM <- ens_2_sym$external_gene_name
  names(gene_syms_BM) <- ens_2_sym$ensembl_gene_id

  test_gene_syms(ens2sym(bbc_obj, BMdataset="mmusculus_gene_ensembl"),
                 gene_syms_BM)
})


###-----------------------------------------------------------------------------

test_that("ens2entrez works", {
  test_entrez_names <- function(bbc, ens2ent){

    bbc_rowdata <- rowData(bbc)

    # test non-na entrez ids correct
    nonna_genes <- bbc_rowdata[!is.na(bbc_rowdata$entrez_ids), "entrez_ids", drop=FALSE]
    expect_equivalent(nonna_genes$entrez_ids, ens2ent[rownames(nonna_genes)])

    # test that all non-NA Entrez IDs are unique
    expect_true(length(unique(nonna_genes)) == length(nonna_genes))

    # test that Entrez IDs with more than one Ensembl match were converted to NAs
    ens2ent_in_bbc <- ens2ent[names(ens2ent) %in% rownames(bbc)]
    dup_ens2ent <- ens2ent_in_bbc[duplicated(ens2ent_in_bbc)]
    ens_genes_w_dup_entrez <- names(ens2ent_in_bbc[ens2ent_in_bbc %in% dup_ens2ent])
    bbc_rowdata_dup <- bbc_rowdata[ens_genes_w_dup_entrez, "entrez_ids"]
    expect_true(all(is.na(bbc_rowdata_dup)))

    # test that NAs are either Entrez IDs with matches to more than one ensembl
    # gene or were absent from the database
    na_genes <- bbc_rowdata[is.na(bbc_rowdata$entrez_ids), "entrez_ids",
                                 drop=FALSE]
    expect_equivalent(nrow(na_genes),
                      sum(rownames(bbc_rowdata) %in% ens_genes_w_dup_entrez |
                            !rownames(bbc_rowdata) %in% names(ens2ent_in_bbc)))

    invisible()
  }

  # OrgDb
  ens_genes <-  AnnotationDbi::keys(org.Mm.eg.db, keytype="ENSEMBL")

  entrez_ids <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                      keys = ens_genes,
                                      column = "ENTREZID",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")

  test_entrez_names(ens2entrez(bbc_obj, org.Mm.eg.db), entrez_ids)

  # BioMart
  ensembl <- biomaRt::useMart("ensembl")
  ensembl <- biomaRt::useDataset("mmusculus_gene_ensembl", mart=ensembl)
  ens_2_entrez <- biomaRt::getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = ensembl)

  # take the first match if there are multiple Entrez matches for a particular Ensembl ID
  ens_2_entrez <- ens_2_entrez[!duplicated(ens_2_entrez$ensembl_gene_id), ]

  entrez_ids_BM <- as.character(ens_2_entrez$entrezgene_id)
  names(entrez_ids_BM) <- ens_2_entrez$ensembl_gene_id

  test_entrez_names(ens2entrez(bbc_obj, BMdataset="mmusculus_gene_ensembl"),
                    entrez_ids_BM)

})

###-----------------------------------------------------------------------------



