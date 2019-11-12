
test_that("ens2sym works", {
  bbc_obj2 <- ens2sym(bbc_obj, org.Mm.eg.db)
  bbc_obj2_rowdata <- rowData(bbc_obj2)
  #bbc_obj2_rowdata_isNA <- is.na(bbc_obj2_rowdata$uniq_syms)
  bbc_obj2_rowdata_uniquified <- stringr::str_detect(bbc_obj2_rowdata$uniq_syms,
                                                     "_ENS[A-Z]*G[0-9]+$")

  ens_genes <-  AnnotationDbi::keys(org.Mm.eg.db, keytype="ENSEMBL")

  gene_syms <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                     keys = ens_genes,
                                     column = "SYMBOL",
                                     keytype = "ENSEMBL",
                                     multiVals = "asNA")


  # add genes that are absent from org.db
  missing_genes <- rownames(bbc_obj2_rowdata)[!rownames(bbc_obj2_rowdata) %in%
                                                names(gene_syms)]
  names(missing_genes) <- missing_genes

  gene_syms_df <- data.frame(symbols = c(gene_syms,
                                         rep(NA, length(missing_genes))),
                             stringsAsFactors = FALSE)

  rownames(gene_syms_df) <- c(names(gene_syms), names(missing_genes))

  # replace NAs with Ensembl ID
  gene_syms_df$symbols <- ifelse(is.na(gene_syms_df$symbols),
                                 rownames(gene_syms_df), gene_syms_df$symbols)

  # concatenate gene symbols with Ensembl IDs
  gene_syms_df$concat_sym <- paste0(gene_syms_df$symbols, "_",
                                    rownames(gene_syms_df))

  # match gene order
  stopifnot(all(rownames(bbc_obj2_rowdata) %in% rownames(gene_syms_df)))
  gene_syms_df <- gene_syms_df[rownames(bbc_obj2_rowdata), ]

  # test non-uniquified symbols correct
  expect_equivalent(bbc_obj2_rowdata[!bbc_obj2_rowdata_uniquified, "uniq_syms"],
                    gene_syms_df[!bbc_obj2_rowdata_uniquified, "symbols"])

  # test uniquified symbols correct
  expect_equivalent(bbc_obj2_rowdata[bbc_obj2_rowdata_uniquified, "uniq_syms"],
                    gene_syms_df[bbc_obj2_rowdata_uniquified, "concat_sym"])

  # test uniquified symbols are duplicated
  expect_true(all(table(gene_syms_df[bbc_obj2_rowdata_uniquified, "symbols"]) > 1))
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



