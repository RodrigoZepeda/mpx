rm(list = ls())
pacman::p_load(dplyr, readr, ggplot2, lubridate, scales, ggtext)

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

mxplot <- ggplot(mpx) +
  geom_col(aes(x = date, y = Casos), fill = "#d74f7c") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x        = "",
    y        = "",
    title    = "Viruela Símica",
    subtitle = "Nuevos casos en México (incidencia semanal)",
    caption  = "<br>**Fuente:** Our World in Data | **Github:** @RodrigoZepeda/mpx"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(color = "#d74f7c", hjust = 1, size = 30, 
                                    margin = margin(t = 10, r = 100, b = 0, l = 0, unit = "pt")),
    plot.background  = element_rect(fill = "#45485c", color = "#45485c"),
    panel.background = element_rect(fill = "#45485c", color = "#45485c"),
    panel.grid       = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.x     = element_line(color = "#d74f7c"),
    axis.text.x      = element_text(color = "#d74f7c"),
    plot.subtitle    = element_markdown(color = "#7e9bed", hjust = 1),
    plot.caption     = element_markdown(color = "#7e9bed", hjust = 0.5)
  ) +
  geom_text(aes(x = date, y = Casos + 1, label = scales::comma(Casos)), vjust = 0,
            size = 3, color = "#d74f7c")
ggsave("images/Monkeypox_Mx.png", plot = mxplot, dpi = 750, width = 6, height = 4)
ggsave("images/Monkeypox_Mx.pdf", plot = mxplot, dpi = 750, width = 6, height = 4)


