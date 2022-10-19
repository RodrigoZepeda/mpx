rm(list = ls())
pacman::p_load(dplyr, readr, ggplot2, cmdstanr, lubridate, stringr)
options(mc.cores = parallel::detectCores() - 2)

mpx <- readr::read_rds("mpx.rds")

# total count
N <- 1.e5;

# times
n_days <- length(mpx$Casos) 

initial_values <- c("S"   = N - mpx$Casos[1], 
                    "I"   = mpx$Casos[1], 
                    "R"   = 0,
                    "PSS" = 0, 
                    "PSI" = 0, 
                    "PII" = 0, 
                    "PSR" = 0, 
                    "PIR" = 0, 
                    "PRR" = 0, 
                    "T"   = mpx$Casos[1])

# data for Stan
data_sir <- list(Nweeks         = length(mpx$Casos), 
                 N              = N, 
                 initial_values = initial_values, 
                 cases          = mpx$Casos)

# number of MCMC steps
niter <- 2000

model <- cmdstanr::cmdstan_model("mpx.stan")

sims  <- model$sample(data = data_sir, iter_warmup = 200, 
                      iter_sampling = 200, chains = 1, seed = 25986)

model_variables <- sims$summary("mpx") |>
  mutate(vartype = stringr::str_remove_all(variable, "mpx\\[[0-9]+,|\\]")) |>
  mutate(vartime = stringr::str_remove_all(variable, "mpx\\[|,[0-9]+\\]")) |>
  mutate(
    vartype = case_when(
      vartype == "1"  ~ "S",
      vartype == "2"  ~ "I",
      vartype == "3"  ~ "R",
      vartype == "4"  ~ "PSS",
      vartype == "5"  ~ "PSI",
      vartype == "6"  ~ "PII",
      vartype == "7"  ~ "PSR",
      vartype == "8"  ~ "PIR",
      vartype == "9"  ~ "PRR",
      vartype == "10" ~ "Total infected",
    )
  ) |>
  mutate(vartime = as.numeric(vartime))

model_variables |>
  #filter(vartype == "Total infected") |>
  ggplot() +
  geom_ribbon(aes(x = vartime, ymin = q5, ymax = q95), fill = "deepskyblue4", alpha = 0.5) +
  geom_line(aes(x = vartime, y = median), color = "deepskyblue4") +
  facet_wrap(~vartype, scales = "free_y") +
  #geom_point(aes(x = vartime, y = Casos), data = mpx |> mutate(vartime = 1:n())) +
  theme_bw()

