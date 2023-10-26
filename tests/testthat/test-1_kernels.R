# CompStatKernel
test_that("Kernel Initializes correctly", {
  g <- CompStatKernel("g")
  e <- CompStatKernel("e")
  expect_equal(class(g), "CompStatKernel")
  expect_equal(class(e), "CompStatKernel")

  g <- CompStatKernel(g)
  e <- CompStatKernel(e)
  expect_equal(class(g), "CompStatKernel")
  expect_equal(class(e), "CompStatKernel")

  g <- CompStatKernel("GGGGGGGaussian")
  e <- CompStatKernel("e______ppaaaanne   ")
  expect_equal(class(g), "CompStatKernel")
  expect_equal(class(e), "CompStatKernel")

  expect_error(CompStatKernel(NULL))
  expect_error(CompStatKernel(1))
})

test_that("Implemented kernels are calculated correctly", {
  grid <- seq(-2, 2, 0.001)
  g <- CompStatKernel("g")
  e <- CompStatKernel("e")

  expect_equal(g$kernel(grid), dnorm(grid))
  expect_equal(e$kernel(grid), ifelse(abs(grid) < 1, 3/4*(1-grid^2), 0))
})

# kernel_density
test_that("Kernel density estimation is consistent with stats::density.
          Be aware that stats::density rescales its kernels.", {
  set.seed(0)
  x <- rnorm(1000)
  bw <- 0.1
  tol <- 0.01
  g <- CompStatKernel("g")
  e <- CompStatKernel("e")
  ref_g <- density(x, bw = bw, kernel = "g")
  ref_e <- density(x, bw = bw, kernel = "e")
  test_g <- kernel_density(g, x, bw, grid = ref_g$x)
  test_e <- kernel_density(e, x, bw / sqrt(e$sigma2), grid = ref_e$x)
  expect_true(all(abs(ref_g$y - test_g) < tol))
  expect_true(all(abs(ref_e$y - test_e) < tol)) # Breaks around (tol < 0.00142)
})

test_that("Kernel density estimates are consistent
          across R and C++ implementations", {
  set.seed(0)
  x <- rnorm(1000)
  bw <- 0.1
  g <- CompStatKernel("g")
  e <- CompStatKernel("e")
  gR <- kernel_density(g, x, bw, method = "r")
  eR <- kernel_density(e, x, bw, method = "r")
  gC <- kernel_density(g, x, bw, method = "c")
  eC <- kernel_density(e, x, bw, method = "c")
  expect_equal(gR, gC)
  expect_equal(eR, eC)
})

# l2norm.CompStatKernel

test_that("L2-Norm estimation works (baseline are matrix implementations)", {
  set.seed(0)
  x <- rnorm(100)
  set.seed(NULL)
  g <- CompStatKernel("g")
  e <- CompStatKernel("e")
  expect_equal(l2norm(g, x, 1), l2norm(g, x, 1, method = "matrix"))
  expect_equal(l2norm(e, x, 1), l2norm(e, x, 1, method = "matrix"))
})
