
# Matrix-method
test_that("Matrix method (singleton)", {
  expect_equal(
    epanechnikov_l2norm_matrix(0, 1),
    4.5
  )
})

test_that("Matrix method (separate singletons)", {
  l <- 8
  x <- seq(4, 4*l, length.out = l)
  expect_equal(
    epanechnikov_l2norm_matrix(x, 1),
    4.5 / l
  )
})

# Can't test the matrix method for 'real' data, since we have nothing to compare
# it to. It serves as our most trusted (and slowest) baseline

# Running
test_that("Running method (singleton)", {
  expect_equal(
    epanechnikov_l2norm_running(0, 1),
    4.5
  )
})

test_that("Running method (separate singletons)", {
  l <- 8
  x <- seq(4, 4*l, length.out = l)
  expect_equal(
    epanechnikov_l2norm_running(x, 1),
    4.5 / l
  )
})

test_that("Running method (real data)", {
  x <- rnorm(1000)
  expect_equal(
    epanechnikov_l2norm_matrix(x, 0.1),
    epanechnikov_l2norm_running(x, 0.1)
  )
})

# Binning
test_that("Binning method (singleton)", {
  expect_equal(
    epanechnikov_l2norm_binning(0, 1),
    4.5
  )
})

test_that("Binning method (separate singletons)", {
  l <- 8
  x <- seq(4, 4*l, length.out = l)
  expect_equal(
    epanechnikov_l2norm_binning(x, 1),
    4.5 / l
  )
})

test_that("Binning method (real data)", {
  x <- rnorm(1000)
  expect_equal(
    epanechnikov_l2norm_matrix(x, 0.1),
    epanechnikov_l2norm_binning(x, 0.1)
  )
})

# Cpp running
test_that("Running cpp method (singleton)", {
  expect_equal(
    epanechnikov_l2norm_runningC(0, 1),
    4.5
  )
})

test_that("Running cpp method (separate singletons)", {
  l <- 8
  x <- seq(4, 4*l, length.out = l)
  expect_equal(
    epanechnikov_l2norm_runningC(x, 1),
    4.5 / l
  )
})

test_that("Running cpp method (real data)", {
  x <- rnorm(1000)
  expect_equal(
    epanechnikov_l2norm_matrix(x, 0.1),
    epanechnikov_l2norm_runningC(x, 0.1)
  )
})



