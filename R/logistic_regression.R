
data(horses)

glm(dead~splines::bs(temp) + natemp, data = horses, family = binomial())

splines::splineDesign(quantile(horses$temp, probs=seq(0.1, 0.9, 0.1)), horses$temp, outer.ok = T)

splines::bs(horses$temp)

design <- matrix(c(1,0,0,2), nrow = 2)
coef <- c(1,1)

x <- seq(-3,3, 0.001)
y <- poly(x, 5)

data <- tibble::tibble(
  x = x
) %>% cbind(y) %>%
  pivot_longer(cols = -c("x"))

ggplot(data, aes(x = x, y = value, color = name)) +
  geom_line()

logistic_loglikelihood <- function(response, design, coef, lambda = 1){
  linear_predictor <- design %*% coef
  exp_lp <- exp(linear_predictor)
  p <- exp_lp / (1 + exp_lp)
  lambda * mean(response * log(p) + (1-response) * log(1-p))
}
