
# Legacy - Median imputation with missingness indicator for temperature
# horses <- readr::read_csv("data-raw/horses.csv") %>%
#   mutate(temp = Temperature %>% replace_na(replace = Temperature %>% median(na.rm = T)),
#          natemp = is.na(Temperature) %>% as.integer()
#          ) %>%
#   select(-c(Temperature))

horses <- readr::read_csv("data-raw/horses.csv") %>%
  tidyr::drop_na()

usethis::use_data(horses)
