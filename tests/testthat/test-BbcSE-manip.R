test_that("BbcSE constructors return valid BbcSE objects.", {
  expect_true(validObject(bbc_obj))
  expect_true(validObject(BbcSE(counts = counts_mat, aln_rates = aln_rates)))
  expect_true(validObject(BbcSE())) # exported
  expect_true(validObject(.BbcSE())) # internal
  expect_error(BbcSE(counts=1))
  expect_error(BbcSE(aln_rates=1))
  expect_error(BbcSE(granges=1))
})


test_that("aln_rates getter works", {
  expect_identical(rownames(aln_rates(bbc_obj)), colnames(bbc_obj))
  expect_identical(
    aln_rates(bbc_obj, withDimnames=FALSE)[, colnames(aln_rates)], aln_rates)
})

test_that("aln_rates setter works", {
  bbc_obj2 <- bbc_obj # make a new copy of the BbcSE object
  old_aln_rates <- aln_rates(bbc_obj2)
  aln_rates(bbc_obj2) <- old_aln_rates * 2

  expect_equivalent(aln_rates(bbc_obj2), old_aln_rates * 2)
  expect_error(aln_rates(bbc_obj2) <- 0)
})
