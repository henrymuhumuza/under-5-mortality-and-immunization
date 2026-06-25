# ============================================================
# 03_merge_panel.R
# Purpose:
# Merge cleaned World Bank and UNICEF datasets into one SSA panel
# Output:
# data/processed/immunization_panel_ssa_2000_2024.csv
# ============================================================

library(dplyr)
library(readr)

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

cleaned_dir <- file.path(base_dir, "data", "cleaned")
processed_dir <- file.path(base_dir, "data", "processed")

# ----------------------------
# Read cleaned datasets
# ----------------------------

u5mr <- read_csv(file.path(cleaned_dir, "u5mr_clean.csv"), show_col_types = FALSE)
dtp3 <- read_csv(file.path(cleaned_dir, "dtp3_clean.csv"), show_col_types = FALSE)
mcv1 <- read_csv(file.path(cleaned_dir, "mcv1_clean.csv"), show_col_types = FALSE)
pol3 <- read_csv(file.path(cleaned_dir, "pol3_clean.csv"), show_col_types = FALSE)

gdp_pc <- read_csv(file.path(cleaned_dir, "gdp_pc_clean.csv"), show_col_types = FALSE)
health_exp_pc <- read_csv(file.path(cleaned_dir, "health_exp_pc_clean.csv"), show_col_types = FALSE)
fertility <- read_csv(file.path(cleaned_dir, "fertility_clean.csv"), show_col_types = FALSE)
female_secondary <- read_csv(file.path(cleaned_dir, "female_secondary_enrollment_clean.csv"), show_col_types = FALSE)
basic_water <- read_csv(file.path(cleaned_dir, "basic_water_clean.csv"), show_col_types = FALSE)
basic_sanitation <- read_csv(file.path(cleaned_dir, "basic_sanitation_clean.csv"), show_col_types = FALSE)
physicians <- read_csv(file.path(cleaned_dir, "physicians_per_1000_clean.csv"), show_col_types = FALSE)

# ----------------------------
# Sub-Saharan Africa ISO3 list
# ----------------------------

ssa_iso3 <- c(
  "AGO","BEN","BWA","BFA","BDI","CPV","CMR","CAF","TCD","COM",
  "COG","COD","CIV","DJI","GNQ","ERI","SWZ","ETH","GAB","GMB",
  "GHA","GIN","GNB","KEN","LSO","LBR","MDG","MWI","MLI","MRT",
  "MUS","MOZ","NAM","NER","NGA","RWA","STP","SEN","SYC","SLE",
  "SOM","ZAF","SSD","SDN","TZA","TGO","UGA","ZMB","ZWE"
)

# ----------------------------
# Merge global panel first
# ----------------------------

master_global <- u5mr %>%
  left_join(dtp3, by = c("iso3", "year")) %>%
  left_join(mcv1, by = c("iso3", "year")) %>%
  left_join(pol3, by = c("iso3", "year")) %>%
  left_join(gdp_pc, by = c("iso3", "year")) %>%
  left_join(health_exp_pc, by = c("iso3", "year")) %>%
  left_join(fertility, by = c("iso3", "year")) %>%
  left_join(female_secondary, by = c("iso3", "year")) %>%
  left_join(basic_water, by = c("iso3", "year")) %>%
  left_join(basic_sanitation, by = c("iso3", "year")) %>%
  left_join(physicians, by = c("iso3", "year")) %>%
  arrange(iso3, year)

# ----------------------------
# Restrict to SSA
# ----------------------------

master_ssa <- master_global %>%
  filter(iso3 %in% ssa_iso3) %>%
  arrange(iso3, year)

# ----------------------------
# Quality checks
# ----------------------------

duplicate_check <- master_ssa %>%
  count(iso3, year) %>%
  filter(n > 1)

missing_summary <- data.frame(
  variable = names(master_ssa),
  missing = sapply(master_ssa, function(x) sum(is.na(x))),
  pct_missing = round(100 * sapply(master_ssa, function(x) mean(is.na(x))), 1)
) %>%
  arrange(desc(pct_missing))

country_count <- master_ssa %>%
  summarise(
    countries = n_distinct(iso3),
    min_year = min(year),
    max_year = max(year),
    rows = n()
  )

print(country_count)
print(duplicate_check)
print(missing_summary)

# ----------------------------
# Save outputs
# ----------------------------

write_csv(
  master_global,
  file.path(processed_dir, "immunization_panel_global_2000_2024.csv")
)

write_csv(
  master_ssa,
  file.path(processed_dir, "immunization_panel_ssa_2000_2024.csv")
)

write_csv(
  missing_summary,
  file.path(processed_dir, "missing_summary_ssa.csv")
)

message("Merge complete.")