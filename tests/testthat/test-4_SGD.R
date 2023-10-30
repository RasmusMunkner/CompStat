test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})




batch_gradient(
  design = matrix(c(1, 0, 0, 2, 1, 4), ncol = 2),
  coef = c(1,2),
  y = c(1.5,3, 2)
  )

res <- SGD_CPP(design, coef, y, 0.001, 500, 1)

plot(seq_along(res$obj), res$obj)




eta <- design %*% res$coef[[500]]
p <- exp(eta) / (1+exp(eta))
diag_dummy <- diag(1, length(eta))
diag(diag_dummy) <- eta / (1+ eta)^2
dp <-  diag_dummy %*% design
g <- -1 / length(y) * sum(y * log(p) + (1-y) * log(1-p))
dg <- - (y / p + (1-y) / (1-p)) / length(y)
grad <- t(dp) %*% dg
grad
