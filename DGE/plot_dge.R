rm(list = ls())
pacman::p_load(lubridate, ggtext, tidyverse)

#Datos de Our World in Data
mpx <- read_csv("DGE/casos_csv/MPX_reporte_dge.csv",
                col_types = cols(
                  Fecha_reporte  = col_date(format = "%Y-%m-%d"),
                  Estado         = col_character(),
                )) |>
  mutate(Estado = if_else(Estado == "Baja California **","Baja California", Estado)) |>
  mutate(Estado = if_else(Estado == "México","Estado de México", Estado)) |>
  mutate(Estado = if_else(Estado == "Queretaro","Querétaro", Estado))

#Relleno de ceros
mpx <- tibble(
  Fecha_reporte = unique(mpx$Fecha_reporte)
  ) |>
  expand_grid(Estado = unique(mpx$Estado)) |>
  left_join(mpx, by = c("Fecha_reporte", "Estado")) |>
  mutate(Casos = if_else(is.na(Casos), 0, Casos))

mpx <- mpx |>
  arrange(Fecha_reporte, Estado) |>
  group_by(Estado) |>
  mutate(Incidencia = Casos - lag(Casos, default = 0)) |>
  ungroup()

mxplot <- ggplot(mpx) +
  geom_col(aes(x = Fecha_reporte, y = Incidencia), fill = "#d74f7c") +
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
  geom_text(aes(x = Fecha_reporte, y = Incidencia, 
                label = Incidencia), vjust = 0,
            size = 1.2, color = "#d74f7c") +
  facet_wrap(~Estado, ncol = 4)
ggsave("images/Monkeypox_Mx_DGE.png", plot = mxplot, dpi = 750, width = 8, height = 10)
ggsave("images/Monkeypox_Mx_DGE.pdf", plot = mxplot, width = 8, height = 10)

mpx |>
  write_excel_csv("data/MPX_DGE.csv")

mpx |>
  group_by(Fecha_reporte) |>
  summarise(Casos = sum(Casos)) |>
  ggplot() +
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
  geom_text(aes(x = Fecha_reporte, y = Casos + 1, label = Casos), vjust = 0,
            size = 3, color = "#d74f7c")
ggsave("images/Monkeypox_Mx.png",  dpi = 750, width = 6, height = 4)
ggsave("images/Monkeypox_Mx.pdf", dpi = 750, width = 6, height = 4)

