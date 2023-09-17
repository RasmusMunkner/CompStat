## code to prepare `poisson` dataset goes here

poisson <- readr::read_csv("data-raw/poisson.csv")
usethis::use_data(poisson, overwrite = TRUE)
