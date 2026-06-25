# ============================================================
# 02_clean_unicef_antigens.R
# Purpose:
# Clean UNICEF antigen coverage Excel files from wide format to long format
# Output:
# Standardized CSV files with columns: iso3, year, antigen variable
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# ----------------------------
# Project paths
# ----------------------------

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

raw_dir <- file.path(base_dir, "data", "raw")
cleaned_dir <- file.path(base_dir, "data", "cleaned")

# ----------------------------
# Cleaning function
# ----------------------------

clean_unicef_antigen <- function(input_file, value_name, output_file) {
  
  message("Cleaning: ", input_file)
  
  raw <- read_excel(input_file)
  
  clean <- raw %>%
    filter(!is.na(iso3)) %>%
    filter(str_detect(iso3, "^[A-Z]{3}$")) %>%
    select(iso3, matches("^\\d{4}$")) %>%
    pivot_longer(
      cols = -iso3,
      names_to = "year",
      values_to = "value"
    ) %>%
    mutate(
      year = as.integer(year),
      value = as.numeric(value)
    ) %>%
    filter(year >= 2000, year <= 2024) %>%
    rename(!!value_name := value) %>%
    arrange(iso3, year)
  
  write_csv(clean, output_file)
  
  message("Saved: ", output_file)
  message("Rows: ", nrow(clean))
  message("")
  
  return(clean)
}

# ----------------------------
# Clean UNICEF antigen datasets
# ----------------------------

dtp3 <- clean_unicef_antigen(
  input_file = file.path(raw_dir, "dtp3.xlsx"),
  value_name = "dtp3",
  output_file = file.path(cleaned_dir, "dtp3_clean.csv")
)

mcv1 <- clean_unicef_antigen(
  input_file = file.path(raw_dir, "mcv1.xlsx"),
  value_name = "mcv1",
  output_file = file.path(cleaned_dir, "mcv1_clean.csv")
)

pol3 <- clean_unicef_antigen(
  input_file = file.path(raw_dir, "pol3.xlsx"),
  value_name = "pol3",
  output_file = file.path(cleaned_dir, "pol3_clean.csv")
)

# ----------------------------
# Quality checks
# ----------------------------

unicef_datasets <- list(
  dtp3 = dtp3,
  mcv1 = mcv1,
  pol3 = pol3
)

quality_check <- lapply(names(unicef_datasets), function(name) {
  
  df <- unicef_datasets[[name]]
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
  file.path(cleaned_dir, "unicef_cleaning_quality_check.csv")
)

message("UNICEF antigen cleaning complete.")