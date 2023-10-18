
data(horses)

glm(dead~splines::bs(Temperature), data = horses, family = binomial())

design <- matrix(c(1,0,0,2), nrow = 2)
coef <- c(1,1)

x <- seq(-3,3, 0.001)
y <- poly(x, 5)

data <- tibble::tibble(
  x = x
) %>% cbind(y) %>%
  tidyr::pivot_longer(cols = -c("x"))

ggplot2::ggplot(data, ggplot2::aes(x = x, y = value, color = name)) +
  ggplot2::geom_line()


