test_that("BbcSE constructors return valid BbcSE objects.", {
  expect_true(validObject(bbc_obj))
  expect_true(validObject(BbcSE(counts = counts_mat, aln_rates = aln_rates)))
  expect_true(validObject(BbcSE())) # exported
  expect_true(validObject(.BbcSE())) # internal
})
