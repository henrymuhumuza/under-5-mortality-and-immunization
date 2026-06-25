# ============================================================
# 07_model_diagnostics_and_sensitivity.R
# Purpose:
# Diagnose model quality and run sensitivity analyses
# Outputs:
# - VIF table
# - Separate antigen models
# - Sensitivity models excluding high-income/small SSA countries
# - Lagged DTP3 model
# - Diagnostic plots
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(fixest)
library(broom)
library(car)
library(modelsummary)

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

processed_dir <- file.path(base_dir, "data", "processed")
tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")

# ----------------------------
# Read data
# ----------------------------

immunization_panel <- read_csv(
  file.path(processed_dir, "immunization_panel_ssa_2000_2024.csv"),
  show_col_types = FALSE
)

# ----------------------------
# Main complete-case dataset
# ----------------------------

main_df <- immunization_panel %>%
  select(
    iso3,
    year,
    u5mr,
    dtp3,
    mcv1,
    pol3,
    gdp_pc,
    health_exp_pc,
    fertility,
    basic_water,
    basic_sanitation
  ) %>%
  drop_na() %>%
  mutate(
    log_gdp_pc = log(gdp_pc),
    log_health_exp_pc = log(health_exp_pc)
  )

# ============================================================
# 1. Multicollinearity check
# ============================================================

vif_model <- lm(
  u5mr ~ dtp3 + mcv1 + pol3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = main_df
)

vif_table <- data.frame(
  variable = names(car::vif(vif_model)),
  vif = as.numeric(car::vif(vif_model))
) %>%
  arrange(desc(vif))

write_csv(
  vif_table,
  file.path(tables_dir, "vif_table_multiple_antigen_model.csv")
)

print(vif_table)

# ============================================================
# 2. Separate antigen models
# ============================================================

model_dtp3 <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = main_df,
  cluster = ~ iso3
)

model_mcv1 <- feols(
  u5mr ~ mcv1 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = main_df,
  cluster = ~ iso3
)

model_pol3 <- feols(
  u5mr ~ pol3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = main_df,
  cluster = ~ iso3
)

separate_antigen_models <- list(
  "DTP3 model" = model_dtp3,
  "MCV1 model" = model_mcv1,
  "Pol3 model" = model_pol3
)

modelsummary(
  separate_antigen_models,
  output = file.path(tables_dir, "separate_antigen_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# 3. Fixed-effects separate antigen models
# ============================================================

model_dtp3_fe <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3 + year,
  data = main_df,
  cluster = ~ iso3
)

model_mcv1_fe <- feols(
  u5mr ~ mcv1 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3 + year,
  data = main_df,
  cluster = ~ iso3
)

model_pol3_fe <- feols(
  u5mr ~ pol3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3 + year,
  data = main_df,
  cluster = ~ iso3
)

separate_antigen_fe_models <- list(
  "DTP3 FE" = model_dtp3_fe,
  "MCV1 FE" = model_mcv1_fe,
  "Pol3 FE" = model_pol3_fe
)

modelsummary(
  separate_antigen_fe_models,
  output = file.path(tables_dir, "separate_antigen_fixed_effects_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# 4. Sensitivity analysis excluding wealthier/small SSA countries
# ============================================================

exclude_iso3 <- c(
  "MUS", # Mauritius
  "SYC", # Seychelles
  "ZAF", # South Africa
  "BWA", # Botswana
  "NAM", # Namibia
  "GAB"  # Gabon
)

restricted_df <- main_df %>%
  filter(!iso3 %in% exclude_iso3)

model_dtp3_restricted <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = restricted_df,
  cluster = ~ iso3
)

model_dtp3_restricted_fe <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3 + year,
  data = restricted_df,
  cluster = ~ iso3
)

sensitivity_models <- list(
  "Main adjusted" = model_dtp3,
  "Main FE" = model_dtp3_fe,
  "Restricted adjusted" = model_dtp3_restricted,
  "Restricted FE" = model_dtp3_restricted_fe
)

modelsummary(
  sensitivity_models,
  output = file.path(tables_dir, "sensitivity_excluding_high_income_ssa.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# 5. Lagged DTP3 model
# Question: Does DTP3 coverage predict next-year U5MR?
# ============================================================

lagged_df <- main_df %>%
  arrange(iso3, year) %>%
  group_by(iso3) %>%
  mutate(
    u5mr_next_year = lead(u5mr, 1)
  ) %>%
  ungroup() %>%
  drop_na(u5mr_next_year)

model_dtp3_lagged <- feols(
  u5mr_next_year ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = lagged_df,
  cluster = ~ iso3
)

model_dtp3_lagged_fe <- feols(
  u5mr_next_year ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3 + year,
  data = lagged_df,
  cluster = ~ iso3
)

lagged_models <- list(
  "Lagged adjusted" = model_dtp3_lagged,
  "Lagged FE" = model_dtp3_lagged_fe
)

modelsummary(
  lagged_models,
  output = file.path(tables_dir, "lagged_dtp3_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# 6. Combine key DTP3 estimates for forest plot
# ============================================================

forest_df <- bind_rows(
  tidy(model_dtp3, conf.int = TRUE) %>% mutate(model = "Adjusted"),
  tidy(model_dtp3_fe, conf.int = TRUE) %>% mutate(model = "Country + year FE"),
  tidy(model_dtp3_restricted, conf.int = TRUE) %>% mutate(model = "Restricted adjusted"),
  tidy(model_dtp3_restricted_fe, conf.int = TRUE) %>% mutate(model = "Restricted FE"),
  tidy(model_dtp3_lagged, conf.int = TRUE) %>% mutate(model = "Lagged adjusted"),
  tidy(model_dtp3_lagged_fe, conf.int = TRUE) %>% mutate(model = "Lagged FE")
) %>%
  filter(term == "dtp3") %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "Adjusted",
        "Country + year FE",
        "Restricted adjusted",
        "Restricted FE",
        "Lagged adjusted",
        "Lagged FE"
      )
    )
  )

write_csv(
  forest_df,
  file.path(tables_dir, "dtp3_sensitivity_forest_coefficients.csv")
)

fig_sensitivity_forest <- ggplot(forest_df, aes(x = estimate, y = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18,
    linewidth = 0.8
  ) +
  labs(
    title = "Sensitivity analysis of the DTP3–under-five mortality association",
    subtitle = "Regression coefficients with 95% confidence intervals",
    x = "Change in U5MR per 1 percentage-point increase in DTP3 coverage",
    y = NULL
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure10_dtp3_sensitivity_forest_plot.png"),
  plot = fig_sensitivity_forest,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# 7. Residual diagnostic plots for main adjusted model
# ============================================================

diagnostic_df <- main_df %>%
  mutate(
    fitted_values = fitted(model_dtp3),
    residuals = resid(model_dtp3)
  )

fig_resid_fitted <- ggplot(diagnostic_df, aes(x = fitted_values, y = residuals)) +
  geom_point(alpha = 0.45) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Residuals versus fitted values",
    subtitle = "Main adjusted DTP3 model",
    x = "Fitted values",
    y = "Residuals"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "diagnostic_residuals_vs_fitted.png"),
  plot = fig_resid_fitted,
  width = 8,
  height = 5,
  dpi = 600
)

fig_resid_hist <- ggplot(diagnostic_df, aes(x = residuals)) +
  geom_histogram(bins = 35) +
  labs(
    title = "Distribution of residuals",
    subtitle = "Main adjusted DTP3 model",
    x = "Residuals",
    y = "Count"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "diagnostic_residual_histogram.png"),
  plot = fig_resid_hist,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# 8. Save model objects for later use
# ============================================================

saveRDS(
  list(
    model_dtp3 = model_dtp3,
    model_mcv1 = model_mcv1,
    model_pol3 = model_pol3,
    model_dtp3_fe = model_dtp3_fe,
    model_mcv1_fe = model_mcv1_fe,
    model_pol3_fe = model_pol3_fe,
    model_dtp3_restricted = model_dtp3_restricted,
    model_dtp3_restricted_fe = model_dtp3_restricted_fe,
    model_dtp3_lagged = model_dtp3_lagged,
    model_dtp3_lagged_fe = model_dtp3_lagged_fe
  ),
  file.path(processed_dir, "regression_models.rds")
)

message("Diagnostics and sensitivity analyses complete.")