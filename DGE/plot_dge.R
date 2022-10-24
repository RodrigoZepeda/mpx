rm(list = ls())
pacman::p_load(dplyr, readr, ggplot2, lubridate, scales, ggtext, tidyverse)

#Datos de Our World in Data
mpx <- read_csv("DGE/casos_csv/MPX_reporte_dge.csv",
                col_types = cols(
                  Fecha_reporte  = col_date(format = "%Y-%m-%d"),
                  Estado         = col_character(),
                )) |>
  mutate(Estado = if_else(Estado == "Baja California **","Baja California", Estado))

#Relleno de ceros
mpx <- tibble(
  Fecha_reporte = unique(mpx$Fecha_reporte)
  ) |>
  expand_grid(Estado = unique(mpx$Estado)) |>
  left_join(mpx) |>
  mutate(Casos = if_else(is.na(Casos), 0, Casos))

mxplot <- ggplot(mpx) +
  geom_col(aes(x = Fecha_reporte, y = Casos), fill = "#d74f7c") +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x        = "",
    y        = "",
    title    = "Viruela Símica",
    subtitle = "Casos acumulados en México",
    caption  = "<br>**Fuente:** DGE | **Github:** @RodrigoZepeda/mpx"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(color = "#d74f7c", hjust = 1, size = 50, 
                                    margin = margin(t = 10, r = 100, b = 0, l = 0, unit = "pt")),
    plot.background  = element_rect(fill = "#45485c", color = "#45485c"),
    panel.background = element_rect(fill = "#45485c", color = "#45485c"),
    panel.grid       = element_blank(),
    axis.text.y      = element_blank(),
    strip.text.x     = element_text(color = "#7e9bed"),
    axis.ticks.x     = element_line(color = "#7e9bed"),
    axis.text.x      = element_text(color = "#7e9bed"),
    plot.subtitle    = element_markdown(color = "#7e9bed", hjust = 1, size = 20),
    plot.caption     = element_markdown(color = "#d74f7c", hjust = 0.5)
  ) +
  geom_text(aes(x = Fecha_reporte, y = Casos, 
                label = Casos), vjust = 0,
            size = 1.2, color = "#d74f7c") +
  facet_wrap(~Estado, ncol = 4)
ggsave("images/Monkeypox_Mx_DGE.png", plot = mxplot, dpi = 750, width = 8, height = 10)
ggsave("images/Monkeypox_Mx_DGE.pdf", plot = mxplot, width = 8, height = 10)


