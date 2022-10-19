rm(list = ls())
pacman::p_load(dplyr, readr, ggplot2, cmdstanr, lubridate)

#Datos de Our World in Data
mpx <- read_csv("https://raw.githubusercontent.com/owid/monkeypox/main/owid-monkeypox-data.csv",
                col_types = cols(
                  date     = col_date(format = "%Y-%m-%d"),
                  location = col_character(),
                  iso_code = col_character(),
                  .default = col_double()
                )) |> 
  dplyr::filter(location == "Mexico") |>
  dplyr::select(date, new_cases, new_cases_per_million) |>
  dplyr::mutate(epiweek = epiweek(date)) |>
  dplyr::mutate(epiyear = epiyear(date)) |>
  dplyr::group_by(epiweek, epiyear) |>
  summarise(Casos = sum(new_cases), 
            Casos_por_millon = sum(new_cases_per_million), 
            .groups = "drop") |>
  dplyr::left_join(
    tibble(date = seq(ymd("2022/01/01"), today(), by = "1 day")) |>
      dplyr::mutate(epiweek = epiweek(date)) |>
      dplyr::mutate(epiyear = epiyear(date)) |>
      dplyr::distinct(epiweek, epiyear, .keep_all = TRUE),
    by = c("epiweek","epiyear")
  ) |> 
  dplyr::arrange(date) |>
  readr::write_rds(file = "mpx.rds")


