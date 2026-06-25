# ============================================================
# 09_final_tables.R
# Purpose:
# Generate final manuscript and supplementary tables
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(broom)
library(fixest)
library(modelsummary)
library(psych)

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

processed_dir <- file.path(base_dir, "data", "processed")
tables_dir <- file.path(base_dir, "tables")

# ----------------------------
# Read data
# ----------------------------

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
# Table 1: Descriptive characteristics
# ============================================================

table1 <- main_df %>%
  select(
    u5mr, dtp3, mcv1, pol3,
    gdp_pc, health_exp_pc, fertility,
    basic_water, basic_sanitation
  ) %>%
  psych::describe() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("variable") %>%
  select(variable, n, mean, sd, median, min, max) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(
  table1,
  file.path(tables_dir, "table1_descriptive_characteristics.csv")
)

# ============================================================
# Table 2: Spearman correlation matrix
# ============================================================

corr_vars <- main_df %>%
  select(
    u5mr, dtp3, mcv1, pol3,
    gdp_pc, health_exp_pc, fertility,
    basic_water, basic_sanitation
  )

corr_matrix <- cor(
  corr_vars,
  use = "pairwise.complete.obs",
  method = "spearman"
)

table2 <- as.data.frame(round(corr_matrix, 2)) %>%
  tibble::rownames_to_column("variable")

write_csv(
  table2,
  file.path(tables_dir, "table2_spearman_correlation_matrix.csv")
)

# ============================================================
# Table 3: Separate adjusted antigen models
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

adjusted_models <- list(
  "DTP3 adjusted" = model_dtp3,
  "MCV1 adjusted" = model_mcv1,
  "Pol3 adjusted" = model_pol3
)

modelsummary(
  adjusted_models,
  output = file.path(tables_dir, "table3_adjusted_antigen_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# Table 4: Country + year fixed-effects antigen models
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

fe_models <- list(
  "DTP3 FE" = model_dtp3_fe,
  "MCV1 FE" = model_mcv1_fe,
  "Pol3 FE" = model_pol3_fe
)

modelsummary(
  fe_models,
  output = file.path(tables_dir, "table4_fixed_effects_antigen_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# Table 5: Sensitivity analyses for DTP3
# ============================================================

exclude_iso3 <- c("MUS", "SYC", "ZAF", "BWA", "NAM", "GAB")

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

lagged_df <- main_df %>%
  arrange(iso3, year) %>%
  group_by(iso3) %>%
  mutate(u5mr_next_year = lead(u5mr, 1)) %>%
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

sensitivity_models <- list(
  "Main adjusted" = model_dtp3,
  "Main FE" = model_dtp3_fe,
  "Restricted adjusted" = model_dtp3_restricted,
  "Restricted FE" = model_dtp3_restricted_fe,
  "Lagged adjusted" = model_dtp3_lagged,
  "Lagged FE" = model_dtp3_lagged_fe
)

modelsummary(
  sensitivity_models,
  output = file.path(tables_dir, "table5_dtp3_sensitivity_models.csv"),
  statistic = "conf.int",
  stars = TRUE,
  gof_omit = "IC|Log|Adj|Within|Pseudo"
)

# ============================================================
# Supplementary Table S1: Missingness
# ============================================================

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
  file.path(tables_dir, "supplementary_table_s1_missingness.csv")
)

# ============================================================
# Supplementary Table S2: VIF values
# ============================================================

vif_path <- file.path(tables_dir, "vif_table_multiple_antigen_model.csv")

if (file.exists(vif_path)) {
  vif_table <- read_csv(vif_path, show_col_types = FALSE)
  
  write_csv(
    vif_table,
    file.path(tables_dir, "supplementary_table_s2_vif_values.csv")
  )
}

# ============================================================
# Supplementary Table S3: Country list
# ============================================================

country_list <- immunization_panel %>%
  distinct(iso3) %>%
  arrange(iso3)

write_csv(
  country_list,
  file.path(tables_dir, "supplementary_table_s3_country_list.csv")
)

# ============================================================
# Supplementary Table S4: Tidy coefficients
# ============================================================

all_models_tidy <- bind_rows(
  tidy(model_dtp3, conf.int = TRUE) %>% mutate(model = "DTP3 adjusted"),
  tidy(model_mcv1, conf.int = TRUE) %>% mutate(model = "MCV1 adjusted"),
  tidy(model_pol3, conf.int = TRUE) %>% mutate(model = "Pol3 adjusted"),
  tidy(model_dtp3_fe, conf.int = TRUE) %>% mutate(model = "DTP3 FE"),
  tidy(model_mcv1_fe, conf.int = TRUE) %>% mutate(model = "MCV1 FE"),
  tidy(model_pol3_fe, conf.int = TRUE) %>% mutate(model = "Pol3 FE"),
  tidy(model_dtp3_restricted, conf.int = TRUE) %>% mutate(model = "DTP3 restricted"),
  tidy(model_dtp3_restricted_fe, conf.int = TRUE) %>% mutate(model = "DTP3 restricted FE"),
  tidy(model_dtp3_lagged, conf.int = TRUE) %>% mutate(model = "DTP3 lagged"),
  tidy(model_dtp3_lagged_fe, conf.int = TRUE) %>% mutate(model = "DTP3 lagged FE")
)

write_csv(
  all_models_tidy,
  file.path(tables_dir, "supplementary_table_s4_all_model_coefficients.csv")
)

# ============================================================
# Save model objects
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
  file.path(processed_dir, "final_model_objects.rds")
)

message("Final tables generated successfully.")