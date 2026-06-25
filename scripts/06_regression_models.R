# ============================================================
# 06_regression_models.R
# Purpose:
# Fit regression models for immunization coverage and U5MR
# Outputs:
# - Regression tables
# - Model coefficient files
# - Forest plot
# ============================================================

library(dplyr)
library(readr)
library(broom)
library(fixest)
library(modelsummary)
library(ggplot2)

# ----------------------------
# Project paths
# ----------------------------

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
  drop_na()

# ----------------------------
# Log-transform skewed economic variables
# ----------------------------

main_df <- main_df %>%
  mutate(
    log_gdp_pc = log(gdp_pc),
    log_health_exp_pc = log(health_exp_pc)
  )

# ----------------------------
# Model 1: Unadjusted DTP3
# ----------------------------

model1 <- feols(
  u5mr ~ dtp3,
  data = main_df,
  cluster = ~ iso3
)

# ----------------------------
# Model 2: Adjusted DTP3
# ----------------------------

model2 <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = main_df,
  cluster = ~ iso3
)

# ----------------------------
# Model 3: Multiple antigen model
# ----------------------------

model3 <- feols(
  u5mr ~ dtp3 + mcv1 + pol3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation,
  data = main_df,
  cluster = ~ iso3
)

# ----------------------------
# Model 4: Country fixed effects
# ----------------------------

model4 <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3,
  data = main_df,
  cluster = ~ iso3
)

# ----------------------------
# Model 5: Country + year fixed effects
# ----------------------------

model5 <- feols(
  u5mr ~ dtp3 + log_gdp_pc + log_health_exp_pc +
    fertility + basic_water + basic_sanitation | iso3 + year,
  data = main_df,
  cluster = ~ iso3
)

# ----------------------------
# Save model summaries
# ----------------------------

models <- list(
  "Unadjusted" = model1,
  "Adjusted" = model2,
  "Multiple antigens" = model3,
  "Country FE" = model4,
  "Country + year FE" = model5
)

modelsummary(
  models,
  output = file.path(tables_dir, "table3_regression_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ----------------------------
# Save tidy coefficients
# ----------------------------

coef_df <- bind_rows(
  tidy(model1, conf.int = TRUE) %>% mutate(model = "Unadjusted"),
  tidy(model2, conf.int = TRUE) %>% mutate(model = "Adjusted"),
  tidy(model3, conf.int = TRUE) %>% mutate(model = "Multiple antigens"),
  tidy(model4, conf.int = TRUE) %>% mutate(model = "Country FE"),
  tidy(model5, conf.int = TRUE) %>% mutate(model = "Country + year FE")
)

write_csv(
  coef_df,
  file.path(tables_dir, "regression_coefficients_tidy.csv")
)

# ----------------------------
# Forest plot: DTP3 coefficient across models
# ----------------------------

dtp3_forest <- coef_df %>%
  filter(term == "dtp3") %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "Unadjusted",
        "Adjusted",
        "Multiple antigens",
        "Country FE",
        "Country + year FE"
      )
    )
  )

fig_forest <- ggplot(
  dtp3_forest,
  aes(x = estimate, y = model)
) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18,
    linewidth = 0.8
  ) +
  labs(
    title = "Association between DTP3 coverage and under-five mortality",
    subtitle = "Regression coefficients with 95% confidence intervals across model specifications",
    x = "Change in U5MR per 1 percentage-point increase in DTP3 coverage",
    y = NULL
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure9_dtp3_forest_plot.png"),
  plot = fig_forest,
  width = 8,
  height = 5,
  dpi = 600
)

# ----------------------------
# Save analysis dataset
# ----------------------------

write_csv(
  main_df,
  file.path(processed_dir, "main_model_complete_case_dataset.csv")
)

# ----------------------------
# Print summaries
# ----------------------------

print(summary(model1))
print(summary(model2))
print(summary(model3))
print(summary(model4))
print(summary(model5))

message("Regression models complete.")