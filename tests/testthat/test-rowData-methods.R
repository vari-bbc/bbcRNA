
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

test_that("ens2entrez works", {
  bbc_obj2 <- ens2entrez(bbc_obj, org.Mm.eg.db)
  bbc_obj2_rowdata <- rowData(bbc_obj2)
  bbc_obj2_rowdata_isna <- rownames(bbc_obj2_rowdata)[is.na(bbc_obj2_rowdata$entrez_ids)]
  bbc_obj2_rowdata_notna <- rownames(bbc_obj2_rowdata)[!is.na(bbc_obj2_rowdata$entrez_ids)]

  # Ensembl IDs in the OrgDb
  ens_genes <-  AnnotationDbi::keys(org.Mm.eg.db, keytype="ENSEMBL")

  # keep only Ensembl IDs present in the BbcSE object
  ens_genes <- ens_genes[ens_genes %in% rownames(bbc_obj)]

  entrez_ids <- AnnotationDbi::mapIds(org.Mm.eg.db,
                                      keys = ens_genes,
                                      column = "ENTREZID",
                                      keytype = "ENSEMBL",
                                      multiVals = "first")

  dup_entrez_ids <- entrez_ids[duplicated(entrez_ids)]
  ens_genes_w_dup_entrez <- names(entrez_ids[entrez_ids %in% dup_entrez_ids])

  # test non-na entrez ids correct
  expect_equivalent(bbc_obj2_rowdata[bbc_obj2_rowdata_notna, "entrez_ids"],
                    entrez_ids[bbc_obj2_rowdata_notna])

  # test that Entrez IDs with more than one Ensembl match were converted to NAs
  expect_true(all(is.na(bbc_obj2_rowdata[ens_genes_w_dup_entrez, "entrez_ids"])))

  # test that all Entrez IDs are unique
  nonna_entrez <- bbc_obj2_rowdata[bbc_obj2_rowdata_notna, "entrez_ids"]
  expect_true(length(unique(nonna_entrez)) == length(nonna_entrez))

  # test that Entrez ID with unique Ensembl ID matches that were NA are absent
  # from org db
  nondup_na_entrez <- bbc_obj2_rowdata[
    !rownames(bbc_obj2_rowdata) %in% ens_genes_w_dup_entrez &
      is.na(bbc_obj2_rowdata$entrez_ids),
    "entrez_ids", drop=FALSE]
  expect_true(all(!rownames(nondup_na_entrez) %in% names(entrez_ids)))
})
