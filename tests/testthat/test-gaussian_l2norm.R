test_that("The integral is correct for a singleton", {
  expect_equal(gaussian_l2norm_matrix(c(0,0,0), 1, debug=F), 3 / 8 / sqrt(pi))
})

test_that("The integral is correct for multiple copies of the same value", {
  expect_equal(gaussian_l2norm_matrix(rep(0,10), 1), 3 / 8 / sqrt(pi))
})


if (F){
  gaussian_l2norm_matrix(c(0,1000), 1)
  estimate_l2norm(c(0, 1000), kernel_code = "g", 1)


  estimate_l2norm(c(0,200), kernel_code = "g", 1)
  gaussian_l2norm_matrix(1000, 1)


  1:10 %>% purrr::map_dbl(.f = function(n) gaussian_l2norm_matrix(rep(0, n), 1, debug=F)) %>%
    log() %>% diff() %>% exp()


  estimate_l2norm(0, kernel_code = "g", 1) / gaussian_l2norm_matrix(0, 1, debug=F)
  estimate_l2norm(c(0,0), kernel_code = "g", 1) / gaussian_l2norm_matrix(c(0,0), 1, debug=F)
  estimate_l2norm(0, kernel_code = "g", 2) / gaussian_l2norm_matrix(0, 2, debug=F)
  estimate_l2norm(c(0,1), kernel_code = "g", 2) / gaussian_l2norm_matrix(c(0,1), 2, debug=F)

  set.seed(0)
  u <- runif(1000)
  estimate_l2norm(u, kernel_code = "g", 2) / gaussian_l2norm_matrix(u, 2, debug=T)
  rm(u)
  set.seed(NULL)
}







