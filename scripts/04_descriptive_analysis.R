# ============================================================
# 04_descriptive_analysis.R
# Purpose:
# Generate descriptive summaries for the SSA immunization panel
# Outputs:
# - Table 1 descriptive summary
# - Missingness summary
# - Correlation matrix
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(psych)
library(corrplot)

# ----------------------------
# Project paths
# ----------------------------

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

processed_dir <- file.path(base_dir, "data", "processed")
tables_dir <- file.path(base_dir, "tables")

# ----------------------------
# Read SSA panel
# ----------------------------

immunization_panel <- read_csv(
  file.path(processed_dir, "immunization_panel_ssa_2000_2024.csv"),
  show_col_types = FALSE
)

# ----------------------------
# Main analysis variables
# ----------------------------

main_vars <- immunization_panel %>%
  select(
    u5mr,
    dtp3,
    mcv1,
    pol3,
    gdp_pc,
    health_exp_pc,
    fertility,
    basic_water,
    basic_sanitation
  )

# ----------------------------
# Missingness summary
# ----------------------------

missing_summary <- data.frame(
  variable = names(immunization_panel),
  missing = sapply(immunization_panel, function(x) sum(is.na(x))),
  pct_missing = round(
    100 * sapply(immunization_panel, function(x) mean(is.na(x))),
    1
  )
) %>%
  arrange(desc(pct_missing))

write_csv(
  missing_summary,
  file.path(tables_dir, "missing_summary_ssa.csv")
)

print(missing_summary)

# ----------------------------
# Descriptive statistics
# ----------------------------

descriptive_summary <- psych::describe(main_vars) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  select(
    variable,
    n,
    mean,
    sd,
    median,
    min,
    max,
    skew,
    kurtosis
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ round(.x, 2)
    )
  )

write_csv(
  descriptive_summary,
  file.path(tables_dir, "table1_descriptive_summary.csv")
)

print(descriptive_summary)

# ----------------------------
# Yearly trend summaries
# ----------------------------

yearly_summary <- immunization_panel %>%
  group_by(year) %>%
  summarise(
    mean_u5mr = mean(u5mr, na.rm = TRUE),
    median_u5mr = median(u5mr, na.rm = TRUE),
    mean_dtp3 = mean(dtp3, na.rm = TRUE),
    mean_mcv1 = mean(mcv1, na.rm = TRUE),
    mean_pol3 = mean(pol3, na.rm = TRUE),
    mean_gdp_pc = mean(gdp_pc, na.rm = TRUE),
    mean_health_exp_pc = mean(health_exp_pc, na.rm = TRUE),
    mean_fertility = mean(fertility, na.rm = TRUE),
    mean_basic_water = mean(basic_water, na.rm = TRUE),
    mean_basic_sanitation = mean(basic_sanitation, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  yearly_summary,
  file.path(tables_dir, "yearly_summary_ssa.csv")
)

print(head(yearly_summary))

# ----------------------------
# Correlation matrix
# ----------------------------

correlation_matrix <- cor(
  main_vars,
  use = "pairwise.complete.obs",
  method = "spearman"
)

correlation_df <- as.data.frame(correlation_matrix) %>%
  tibble::rownames_to_column("variable")

write_csv(
  correlation_df,
  file.path(tables_dir, "table2_spearman_correlation_matrix.csv")
)

print(round(correlation_matrix, 2))

# ----------------------------
# Complete-case count for main model
# ----------------------------

main_model_df <- immunization_panel %>%
  select(
    iso3,
    year,
    u5mr,
    dtp3,
    gdp_pc,
    health_exp_pc,
    fertility,
    basic_water,
    basic_sanitation
  ) %>%
  drop_na()

complete_case_summary <- tibble(
  total_rows = nrow(immunization_panel),
  complete_case_rows = nrow(main_model_df),
  countries = n_distinct(main_model_df$iso3),
  min_year = min(main_model_df$year),
  max_year = max(main_model_df$year)
)

write_csv(
  complete_case_summary,
  file.path(tables_dir, "main_model_complete_case_summary.csv")
)

print(complete_case_summary)

message("Descriptive analysis complete.")