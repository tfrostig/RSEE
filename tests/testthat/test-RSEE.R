set.seed(999)
example.data <- data.frame('x' = rnorm(1000), "y" = rnorm(1000))

test_that("LOWESS returns right colnames", {
  expect_true(all(colnames(RSEE:::LOWESS(example.data$x, d = 2)) == c('loc', 'v', 'a')))
})


test_that("LOWESS returns right number of columns", {
  expect_equal(ncol(RSEE:::LOWESS(example.data$x, d = 3)), 4)
})


test_that("LOWESS returns warning", {
  expect_warning(RSEE:::LOWESS(example.data$x[1:5], d = 2))
})


test_that("RSEE returns expected values", {
  expect_equal(class(RSEE::RSEEsmooth(example.data)), 'list')
})
