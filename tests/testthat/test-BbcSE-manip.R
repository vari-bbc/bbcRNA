test_that("BbcSE constructors return valid BbcSE objects.", {
  expect_true(validObject(bbc_obj))
  expect_true(validObject(BbcSE(counts = counts_mat,
                                aln_metrics = aln_metrics)))
  expect_true(validObject(BbcSE())) # exported
  expect_true(validObject(.BbcSE())) # internal
  expect_error(BbcSE(counts=1))
  expect_error(BbcSE(aln_metrics=1))
  expect_error(BbcSE(granges=1))
})


test_that("aln_metrics getter works", {
  expect_identical(rownames(aln_metrics(bbc_obj)), colnames(bbc_obj))
  expect_identical(
    aln_metrics(bbc_obj, withDimnames=FALSE)[, colnames(aln_metrics)],
    aln_metrics)
})

test_that("aln_metrics setter works", {
  bbc_obj2 <- bbc_obj # make a new copy of the BbcSE object
  old_aln_metrics <- aln_metrics(bbc_obj2)
  aln_metrics(bbc_obj2) <- old_aln_metrics * 2

  expect_equivalent(aln_metrics(bbc_obj2), old_aln_metrics * 2)
  expect_error(aln_metrics(bbc_obj2) <- 0)
})
