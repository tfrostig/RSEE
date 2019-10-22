example.data <- rbind(c(0, 0, 0),
                      c(0, 1, 0),
                      c(0, 0, 0))

rownames(example.data) <- c(-1, 0, 1)
colnames(example.data) <- c(-1, 0, 1)

test_that("Local maxima are found", {
  expect_equal(drop(RSEE:::FindMax(example.data, 0)), c('x' = 0, 'y' = 0))
})
