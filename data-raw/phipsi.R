## code to prepare `phipsi` dataset goes here

library(readr)
library(dplyr)

phipsi <- read_csv("data-raw/phipsi.csv") %>%
  dplyr::select(phi, psi)

usethis::use_data(phipsi, overwrite = TRUE)
