# ============================================================
# 01_clean_worldbank_files.R
# Purpose:
# Clean World Bank-style Excel files from wide format to long format
# Output:
# Standardized CSV files with columns: iso3, year, variable
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(readr)

# ----------------------------
# Project paths
# ----------------------------

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

raw_dir <- file.path(base_dir, "data", "raw")
cleaned_dir <- file.path(base_dir, "data", "cleaned")

# ----------------------------
# Cleaning function
# ----------------------------

clean_worldbank_file <- function(input_file, value_name, output_file) {
  
  message("Cleaning: ", input_file)
  
  raw <- read_excel(input_file)
  
  clean <- raw %>%
    filter(str_detect(`Country Code`, "^[A-Z]{3}$")) %>%
    select(`Country Code`, matches("^\\d{4} \\[YR\\d{4}\\]$")) %>%
    pivot_longer(
      cols = -`Country Code`,
      names_to = "year",
      values_to = "value"
    ) %>%
    mutate(
      year = as.integer(str_extract(year, "\\d{4}")),
      value = na_if(as.character(value), ".."),
      value = as.numeric(value)
    ) %>%
    filter(year >= 2000, year <= 2024) %>%
    rename(iso3 = `Country Code`) %>%
    rename(!!value_name := value) %>%
    arrange(iso3, year)
  
  write_csv(clean, output_file)
  
  message("Saved: ", output_file)
  message("Rows: ", nrow(clean))
  message("")
  
  return(clean)
}

# ----------------------------
# Clean World Bank datasets
# ----------------------------

u5mr <- clean_worldbank_file(
  input_file = file.path(raw_dir, "u5mr.xlsx"),
  value_name = "u5mr",
  output_file = file.path(cleaned_dir, "u5mr_clean.csv")
)

gdp_pc <- clean_worldbank_file(
  input_file = file.path(raw_dir, "gdp_pc.xlsx"),
  value_name = "gdp_pc",
  output_file = file.path(cleaned_dir, "gdp_pc_clean.csv")
)

health_exp_pc <- clean_worldbank_file(
  input_file = file.path(raw_dir, "health_exp_pc.xlsx"),
  value_name = "health_exp_pc",
  output_file = file.path(cleaned_dir, "health_exp_pc_clean.csv")
)

fertility <- clean_worldbank_file(
  input_file = file.path(raw_dir, "fertility.xlsx"),
  value_name = "fertility",
  output_file = file.path(cleaned_dir, "fertility_clean.csv")
)

female_secondary_enrollment <- clean_worldbank_file(
  input_file = file.path(raw_dir, "female_secondary_enrollment.xlsx"),
  value_name = "female_secondary_enrollment",
  output_file = file.path(cleaned_dir, "female_secondary_enrollment_clean.csv")
)

basic_water <- clean_worldbank_file(
  input_file = file.path(raw_dir, "basic_water.xlsx"),
  value_name = "basic_water",
  output_file = file.path(cleaned_dir, "basic_water_clean.csv")
)

basic_sanitation <- clean_worldbank_file(
  input_file = file.path(raw_dir, "basic_sanitation.xlsx"),
  value_name = "basic_sanitation",
  output_file = file.path(cleaned_dir, "basic_sanitation_clean.csv")
)

physicians_per_1000 <- clean_worldbank_file(
  input_file = file.path(raw_dir, "physicians_per_1000.xlsx"),
  value_name = "physicians_per_1000",
  output_file = file.path(cleaned_dir, "physicians_per_1000_clean.csv")
)

# ----------------------------
# Quality checks
# ----------------------------

worldbank_datasets <- list(
  u5mr = u5mr,
  gdp_pc = gdp_pc,
  health_exp_pc = health_exp_pc,
  fertility = fertility,
  female_secondary_enrollment = female_secondary_enrollment,
  basic_water = basic_water,
  basic_sanitation = basic_sanitation,
  physicians_per_1000 = physicians_per_1000
)

quality_check <- lapply(names(worldbank_datasets), function(name) {
  
  df <- worldbank_datasets[[name]]
  value_col <- setdiff(names(df), c("iso3", "year"))
  
  tibble(
    dataset = name,
    rows = nrow(df),
    countries = n_distinct(df$iso3),
    min_year = min(df$year, na.rm = TRUE),
    max_year = max(df$year, na.rm = TRUE),
    missing_iso3 = sum(is.na(df$iso3)),
    duplicate_keys = sum(duplicated(paste(df$iso3, df$year))),
    missing_values = sum(is.na(df[[value_col]])),
    pct_missing = round(100 * mean(is.na(df[[value_col]])), 1)
  )
}) %>%
  bind_rows()

print(quality_check)

write_csv(
  quality_check,
  file.path(cleaned_dir, "worldbank_cleaning_quality_check.csv")
)

message("World Bank cleaning complete.")