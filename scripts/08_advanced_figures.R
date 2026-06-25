# ============================================================
# 08_advanced_figures.R
# Purpose:
# Generate advanced publication-style figures
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(broom)
library(fixest)
library(scales)
library(grid)

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

processed_dir <- file.path(base_dir, "data", "processed")
tables_dir <- file.path(base_dir, "tables")
figures_dir <- file.path(base_dir, "figures")

immunization_panel <- read_csv(
  file.path(processed_dir, "immunization_panel_ssa_2000_2024.csv"),
  show_col_types = FALSE
)

main_df <- immunization_panel %>%
  select(
    iso3, year, u5mr, dtp3, mcv1, pol3,
    gdp_pc, health_exp_pc, fertility,
    basic_water, basic_sanitation
  ) %>%
  drop_na() %>%
  mutate(
    log_gdp_pc = log(gdp_pc),
    log_health_exp_pc = log(health_exp_pc)
  )

# ============================================================
# 1. Separate antigen adjusted models
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

# ============================================================
# 2. Adjusted partial-effect plots: antigen models
# ============================================================

antigen_forest <- bind_rows(
  tidy(model_dtp3, conf.int = TRUE) %>% mutate(model = "DTP3"),
  tidy(model_mcv1, conf.int = TRUE) %>% mutate(model = "MCV1"),
  tidy(model_pol3, conf.int = TRUE) %>% mutate(model = "Pol3")
) %>%
  filter(term %in% c("dtp3", "mcv1", "pol3")) %>%
  mutate(
    antigen = dplyr::recode(
      term,
      dtp3 = "DTP3",
      mcv1 = "MCV1",
      pol3 = "Pol3"
    )
  )

write_csv(
  antigen_forest,
  file.path(tables_dir, "antigen_forest_coefficients.csv")
)

make_partial_effect <- function(model, exposure, label, color) {
  covariate_means <- main_df %>%
    summarise(
      log_gdp_pc = mean(log_gdp_pc, na.rm = TRUE),
      log_health_exp_pc = mean(log_health_exp_pc, na.rm = TRUE),
      fertility = mean(fertility, na.rm = TRUE),
      basic_water = mean(basic_water, na.rm = TRUE),
      basic_sanitation = mean(basic_sanitation, na.rm = TRUE)
    )

  grid <- tibble(
    coverage = seq(
      min(main_df[[exposure]], na.rm = TRUE),
      max(main_df[[exposure]], na.rm = TRUE),
      length.out = 100
    )
  ) %>%
    mutate(
      !!exposure := coverage,
      log_gdp_pc = covariate_means$log_gdp_pc,
      log_health_exp_pc = covariate_means$log_health_exp_pc,
      fertility = covariate_means$fertility,
      basic_water = covariate_means$basic_water,
      basic_sanitation = covariate_means$basic_sanitation
    )

  model_formula <- as.formula(
    paste0(
      "~ ", exposure, " + log_gdp_pc + log_health_exp_pc + ",
      "fertility + basic_water + basic_sanitation"
    )
  )
  x_matrix <- model.matrix(model_formula, data = grid)
  beta <- coef(model)
  vcov_matrix <- vcov(model)
  x_matrix <- x_matrix[, names(beta), drop = FALSE]

  fit <- as.vector(x_matrix %*% beta)
  se <- sqrt(diag(x_matrix %*% vcov_matrix[names(beta), names(beta)] %*% t(x_matrix)))

  grid %>%
    transmute(
      antigen = label,
      coverage = coverage,
      predicted_u5mr = fit,
      conf.low = fit - 1.96 * se,
      conf.high = fit + 1.96 * se,
      color = color
    )
}

partial_effects <- bind_rows(
  make_partial_effect(model_dtp3, "dtp3", "DTP3", "#1b9e77"),
  make_partial_effect(model_mcv1, "mcv1", "MCV1", "#d95f02"),
  make_partial_effect(model_pol3, "pol3", "Pol3", "#7570b3")
)

fig_antigen_partial_effects <- ggplot(
  partial_effects,
  aes(x = coverage, y = predicted_u5mr, color = antigen, fill = antigen)
) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1.1) +
  facet_wrap(~ antigen, nrow = 1) +
  scale_color_manual(values = c("DTP3" = "#1b9e77", "MCV1" = "#d95f02", "Pol3" = "#7570b3")) +
  scale_fill_manual(values = c("DTP3" = "#1b9e77", "MCV1" = "#d95f02", "Pol3" = "#7570b3")) +
  labs(
    title = "Adjusted predicted under-five mortality by vaccine coverage",
    subtitle = "Separate antigen models; other covariates held at complete-case means",
    x = "Coverage (%)",
    y = "Predicted U5MR (deaths per 1,000 live births)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(
  file.path(figures_dir, "figure11_antigen_partial_effects.png"),
  fig_antigen_partial_effects,
  width = 10,
  height = 5,
  dpi = 600,
  bg = "white"
)

# ============================================================
# 3. Full coefficient forest plot for main DTP3 model
# ============================================================

coef_forest <- tidy(model_dtp3, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_clean = dplyr::recode(
      term,
      dtp3 = "DTP3 coverage",
      log_gdp_pc = "Log GDP per capita",
      log_health_exp_pc = "Log health expenditure per capita",
      fertility = "Total fertility rate",
      basic_water = "Basic water access",
      basic_sanitation = "Basic sanitation access"
    ),
    term_clean = factor(term_clean, levels = rev(unique(term_clean)))
  )

write_csv(
  coef_forest,
  file.path(tables_dir, "main_dtp3_model_coefficients.csv")
)

fig_coef_forest <- ggplot(coef_forest, aes(x = estimate, y = term_clean)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  geom_errorbarh(
    aes(xmin = conf.low, xmax = conf.high),
    height = 0.18,
    linewidth = 0.8
  ) +
  labs(
    title = "Main adjusted model coefficients",
    subtitle = "Outcome: under-five mortality rate per 1,000 live births",
    x = "Regression coefficient",
    y = NULL
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(figures_dir, "figure12_main_model_coefficient_forest_plot.png"),
  fig_coef_forest,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# 4. Predicted versus observed U5MR
# ============================================================

prediction_df <- main_df %>%
  mutate(
    predicted_u5mr = fitted(model_dtp3)
  )

fig_pred_obs <- ggplot(prediction_df, aes(x = predicted_u5mr, y = u5mr)) +
  geom_point(alpha = 0.45) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = "Predicted versus observed under-five mortality",
    subtitle = "Main adjusted DTP3 model",
    x = "Predicted U5MR",
    y = "Observed U5MR"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(figures_dir, "figure13_predicted_vs_observed_u5mr.png"),
  fig_pred_obs,
  width = 7,
  height = 6,
  dpi = 600
)

# ============================================================
# 5. Bubble plot: DTP3, U5MR and GDP
# ============================================================

fig_bubble <- ggplot(main_df, aes(x = dtp3, y = u5mr)) +
  geom_point(aes(size = gdp_pc), alpha = 0.35) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_size_continuous(
    labels = comma,
    range = c(1, 8)
  ) +
  labs(
    title = "DTP3 coverage, under-five mortality and GDP per capita",
    subtitle = "Bubble size represents GDP per capita",
    x = "DTP3 coverage (%)",
    y = "Under-five mortality rate per 1,000 live births",
    size = "GDP per capita"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(figures_dir, "figure14_dtp3_u5mr_gdp_bubble_plot.png"),
  fig_bubble,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# 6. Country trajectories for selected countries
# ============================================================

selected_countries <- c(
  "UGA", "KEN", "TZA", "RWA", "ETH", "NGA", "ZAF", "GHA"
)

trajectory_df <- immunization_panel %>%
  filter(iso3 %in% selected_countries) %>%
  arrange(iso3, year)

fig_country_trajectories <- ggplot(
  trajectory_df,
  aes(x = dtp3, y = u5mr, group = iso3)
) +
  geom_path(
    arrow = arrow(length = unit(0.12, "inches")),
    alpha = 0.8
  ) +
  geom_point(size = 1.7, alpha = 0.8) +
  facet_wrap(~ iso3) +
  labs(
    title = "Country trajectories in DTP3 coverage and under-five mortality",
    subtitle = "Movement from 2000 to 2024 among selected sub-Saharan African countries",
    x = "DTP3 coverage (%)",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  file.path(figures_dir, "figure15_selected_country_trajectories.png"),
  fig_country_trajectories,
  width = 10,
  height = 7,
  dpi = 600
)

# ============================================================
# 7. Violin plot: U5MR by period
# ============================================================

period_df <- immunization_panel %>%
  mutate(
    period = case_when(
      year >= 2000 & year <= 2009 ~ "2000–2009",
      year >= 2010 & year <= 2019 ~ "2010–2019",
      year >= 2020 & year <= 2024 ~ "2020–2024"
    ),
    period = factor(
      period,
      levels = c("2000–2009", "2010–2019", "2020–2024")
    )
  )

fig_violin <- ggplot(period_df, aes(x = period, y = u5mr)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.alpha = 0.4) +
  labs(
    title = "Distribution of under-five mortality by period",
    subtitle = "Violin plots with embedded boxplots",
    x = "Period",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  file.path(figures_dir, "figure16_u5mr_violin_by_period.png"),
  fig_violin,
  width = 8,
  height = 5,
  dpi = 600
)

message("Advanced figures generated successfully.")
