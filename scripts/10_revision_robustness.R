# ============================================================
# 10_revision_robustness.R -- analyses added for major revision
# ============================================================
library(dplyr)
library(readr)
library(tidyr)
library(fixest)
library(broom)

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"
processed_dir <- file.path(base_dir, "data", "processed")
tables_dir <- file.path(base_dir, "tables")
panel <- read_csv(file.path(processed_dir, "immunization_panel_ssa_2000_2024.csv"), show_col_types = FALSE)

needed <- c("u5mr", "dtp3", "gdp_pc", "health_exp_pc", "fertility", "basic_water", "basic_sanitation")
analysis_df <- panel |>
  select(iso3, year, all_of(needed)) |>
  mutate(log_gdp_pc = log(gdp_pc), log_health_exp_pc = log(health_exp_pc))

# Included-versus-excluded observations for the complete-case primary analysis.
sample_comparison <- analysis_df |>
  mutate(sample = if_else(if_all(all_of(needed), ~ !is.na(.x)), "Included in primary model", "Excluded from primary model")) |>
  group_by(sample) |>
  summarise(observations = n(), countries = n_distinct(iso3), min_year = min(year), max_year = max(year),
            across(c(u5mr, dtp3, fertility, gdp_pc, health_exp_pc, basic_water, basic_sanitation),
                   list(mean = ~ mean(.x, na.rm = TRUE), missing = ~ sum(is.na(.x))),
                   .names = "{.col}_{.fn}"), .groups = "drop")
write_csv(sample_comparison, file.path(tables_dir, "supplementary_table_s5_included_excluded_comparison.csv"))

main_df <- analysis_df |> drop_na() |> mutate(year_centered = year - min(year))
main_formula <- u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc + fertility + basic_water + basic_sanitation | iso3 + year

# Country clustering permits arbitrary serial dependence within countries. Driscoll--Kraay
# SEs permit general serial and cross-sectional dependence.
mod_cluster <- feols(main_formula, data = main_df, vcov = ~iso3)
mod_dk <- feols(main_formula, data = main_df, vcov = "DK", panel.id = ~iso3 + year)
mod_two_cluster <- feols(main_formula, data = main_df, vcov = ~iso3 + year)
mod_trends <- feols(u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc + fertility + basic_water + basic_sanitation |
                      iso3 + year + iso3[year_centered], data = main_df, vcov = "DK", panel.id = ~iso3 + year)

# Omit health expenditure to retain observations from 2024, when it was unavailable.
no_hexp_df <- panel |> select(iso3, year, u5mr, dtp3, gdp_pc, fertility, basic_water, basic_sanitation) |>
  drop_na() |> mutate(log_gdp_pc = log(gdp_pc))
mod_no_hexp_matched <- feols(u5mr ~ dtp3 + log_gdp_pc + fertility + basic_water + basic_sanitation | iso3 + year,
                             data = filter(no_hexp_df, year <= 2023), vcov = "DK", panel.id = ~iso3 + year)
mod_no_hexp_all <- feols(u5mr ~ dtp3 + log_gdp_pc + fertility + basic_water + basic_sanitation | iso3 + year,
                         data = no_hexp_df, vcov = "DK", panel.id = ~iso3 + year)

dtp3_row <- function(fit, label) {
  obs <- fit$nobs
  tidy(fit, conf.int = TRUE) |> filter(term == "dtp3") |>
    transmute(model = label, estimate, std.error, conf.low, conf.high, p.value, observations = obs)
}

robustness <- bind_rows(
  dtp3_row(mod_cluster, "Country-clustered SEs (primary two-way FE)"), dtp3_row(mod_dk, "Driscoll-Kraay SEs"),
  dtp3_row(mod_two_cluster, "Two-way clustered SEs"), dtp3_row(mod_trends, "Country-specific linear trends; Driscoll-Kraay SEs"),
  dtp3_row(mod_no_hexp_matched, "Omit health expenditure, 2000-2023"), dtp3_row(mod_no_hexp_all, "Omit health expenditure, 2000-2024"))
write_csv(robustness, file.path(tables_dir, "supplementary_table_s6_panel_robustness.csv"))

# Prospective alternative lag structures; covariates are measured in the exposure year.
lag_results <- lapply(1:3, function(k) {
  dat <- main_df |> arrange(iso3, year) |> group_by(iso3) |> mutate(u5mr_future = lead(u5mr, k)) |> ungroup() |> drop_na(u5mr_future)
  mod <- feols(u5mr_future ~ dtp3 + log_gdp_pc + log_health_exp_pc + fertility + basic_water + basic_sanitation | iso3 + year,
               data = dat, vcov = "DK", panel.id = ~iso3 + year)
  dtp3_row(mod, paste0(k, "-year prospective lag; Driscoll-Kraay SEs")) |> mutate(lag_years = k)
}) |> bind_rows()
write_csv(lag_results, file.path(tables_dir, "supplementary_table_s7_alternative_lags.csv"))

# Test whether a quadratic term materially improves the linear DTP3 specification.
mod_quadratic <- feols(u5mr ~ dtp3 + I(dtp3^2) + log_gdp_pc + log_health_exp_pc + fertility + basic_water + basic_sanitation | iso3 + year,
                       data = main_df, vcov = "DK", panel.id = ~iso3 + year)
quadratic_obs <- mod_quadratic$nobs
linearity <- tidy(mod_quadratic, conf.int = TRUE) |> filter(grepl("dtp3\\^2", term)) |>
  mutate(test = "Quadratic DTP3 term in two-way fixed-effects model", observations = quadratic_obs)
write_csv(linearity, file.path(tables_dir, "supplementary_table_s8_linearity_check.csv"))
message("Revision robustness analyses generated successfully.")
