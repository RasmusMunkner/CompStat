

logistic_loglikelihood <- function(coef,
                                   design = horses$Temperature %>% rescale_covariate() %>% poly(degree = 3),
                                   response = horses$dead,
                                   lambda = 1
                                   ){

  linear_predictor <- design %*% coef
  exp_lp <- exp(linear_predictor)
  p <- exp_lp / (1 + exp_lp)
  -mean(response * log(p) + (1-response) * log(1-p)) + lambda * panalty
}

horses$Temperature %>% rescale_covariate() %>%
  cbind(poly(.,6) %>% as.data.frame()) %>%
  `names<-`(c("Temperature", paste0("Degree", 1:6))) %>%
  pivot_longer(cols = -Temperature, names_to="Basis", values_to="splineval") %>%
  ggplot(aes(x = Temperature, y = splineval, color = Basis)) +
  geom_line()

horses$Temperature %>% rescale_covariate() %>% hist()

horses$Temperature %>% rescale_covariate() %>% poly(10) %>% cov()

htemp <- horses$Temperature %>% rescale_covariate()
spline_htemp <- htemp %>% splines::splineDesign(knots = seq(0.1, 0.9, length.out = 8), x = ., outer.ok = T)


tibble(htemp = htemp) %>%
  cbind(spline_htemp %>% setNames(paste0("S", 1:4))) %>%
  as_tibble() %>%
  pivot_longer(cols = -c(htemp)) %>%
  ggplot(aes(x = htemp, y = value, color = name)) +
  geom_line()


