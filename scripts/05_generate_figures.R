# ============================================================
# 05_generate_figures.R
# Purpose:
# Generate exploratory and manuscript-ready figures
# Outputs:
# PNG figures saved in /figures
# ============================================================

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(scales)
library(corrplot)

# ----------------------------
# Project paths
# ----------------------------

base_dir <- "C:/Users/HP/Desktop/IMMUNIZATION STUDY"

processed_dir <- file.path(base_dir, "data", "processed")
figures_dir <- file.path(base_dir, "figures")
tables_dir <- file.path(base_dir, "tables")

# ----------------------------
# Read data
# ----------------------------

immunization_panel <- read_csv(
  file.path(processed_dir, "immunization_panel_ssa_2000_2024.csv"),
  show_col_types = FALSE
)

# ----------------------------
# Create yearly summary
# ----------------------------

yearly_summary <- immunization_panel %>%
  group_by(year) %>%
  summarise(
    mean_u5mr = mean(u5mr, na.rm = TRUE),
    median_u5mr = median(u5mr, na.rm = TRUE),
    mean_dtp3 = mean(dtp3, na.rm = TRUE),
    mean_mcv1 = mean(mcv1, na.rm = TRUE),
    mean_pol3 = mean(pol3, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  yearly_summary,
  file.path(tables_dir, "yearly_summary_for_figures.csv")
)

# ============================================================
# Figure 1: Mean under-five mortality trend
# ============================================================

fig1 <- ggplot(yearly_summary, aes(x = year, y = mean_u5mr)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  labs(
    title = "Mean under-five mortality rate in sub-Saharan Africa, 2000–2024",
    subtitle = "Country-year ecological panel based on World Bank estimates",
    x = "Five-year period",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure1_mean_u5mr_trend.png"),
  plot = fig1,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# Figure 2: Mean vaccine coverage trends
# ============================================================

vaccine_trends <- yearly_summary %>%
  select(year, mean_dtp3, mean_mcv1, mean_pol3) %>%
  pivot_longer(
    cols = starts_with("mean_"),
    names_to = "antigen",
    values_to = "coverage"
  ) %>%
  mutate(
    antigen = recode(
      antigen,
      mean_dtp3 = "DTP3",
      mean_mcv1 = "MCV1",
      mean_pol3 = "Pol3"
    )
  )

fig2 <- ggplot(vaccine_trends, aes(x = year, y = coverage, color = antigen, linetype = antigen)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  labs(
    title = "Mean routine childhood immunization coverage in sub-Saharan Africa, 2000–2024",
    subtitle = "DTP3, MCV1 and Pol3 coverage estimates from UNICEF",
    x = "Year",
    y = "Coverage (%)",
    color = "Antigen",
    linetype = "Antigen"
  ) +
  ggtitle("Mean routine childhood immunization coverage, 2000-2024",
          subtitle = "DTP3, MCV1 and Pol3 coverage estimates from UNICEF") +
  scale_color_manual(values = c("DTP3" = "#1b9e77", "MCV1" = "#d95f02", "Pol3" = "#7570b3")) +
  scale_y_continuous(limits = c(55, 85), breaks = seq(55, 85, 5)) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure2_vaccine_coverage_trends.png"),
  plot = fig2,
  width = 8,
  height = 5,
  dpi = 600,
  bg = "white"
)

# ============================================================
# Figure 3: Boxplot of U5MR by year
# ============================================================

u5mr_by_period <- immunization_panel %>%
  mutate(
    period = cut(
      year,
      breaks = c(1999, 2004, 2009, 2014, 2019, 2024),
      labels = c("2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2024")
    )
  )

fig3 <- ggplot(u5mr_by_period, aes(x = period, y = u5mr)) +
  geom_boxplot(outlier.alpha = 0.4) +
  labs(
    title = "Distribution of under-five mortality across sub-Saharan African countries",
    subtitle = "Boxplots by year, 2000–2024",
    x = "Year",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 12) +
  labs(subtitle = "Country-year observations grouped into five-year periods, 2000-2024")

ggsave(
  filename = file.path(figures_dir, "figure3_u5mr_boxplot_by_year.png"),
  plot = fig3,
  width = 11,
  height = 6,
  dpi = 600
)

# ============================================================
# Figure 4: DTP3 vs U5MR scatterplot
# ============================================================

fig4 <- ggplot(immunization_panel, aes(x = dtp3, y = u5mr)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Association between DTP3 coverage and under-five mortality",
    subtitle = "Each point represents one country-year observation",
    x = "DTP3 coverage (%)",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure4_dtp3_u5mr_scatter.png"),
  plot = fig4,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# Figure 5: Faceted vaccine boxplots by decade
# ============================================================

vaccine_long <- immunization_panel %>%
  mutate(
    period = case_when(
      year >= 2000 & year <= 2009 ~ "2000–2009",
      year >= 2010 & year <= 2019 ~ "2010–2019",
      year >= 2020 & year <= 2024 ~ "2020–2024"
    )
  ) %>%
  select(year, period, dtp3, mcv1, pol3) %>%
  pivot_longer(
    cols = c(dtp3, mcv1, pol3),
    names_to = "antigen",
    values_to = "coverage"
  ) %>%
  mutate(
    antigen = recode(
      antigen,
      dtp3 = "DTP3",
      mcv1 = "MCV1",
      pol3 = "Pol3"
    )
  )

fig5 <- ggplot(vaccine_long, aes(x = period, y = coverage)) +
  geom_boxplot(outlier.alpha = 0.4) +
  facet_wrap(~ antigen) +
  labs(
    title = "Distribution of vaccine coverage by period",
    subtitle = "Country-year observations grouped into three periods",
    x = "Period",
    y = "Coverage (%)"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure5_vaccine_coverage_boxplots_by_period.png"),
  plot = fig5,
  width = 9,
  height = 5,
  dpi = 600
)

# ============================================================
# Figure 6: Correlation heatmap
# ============================================================

corr_vars <- immunization_panel %>%
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

corr_matrix <- cor(
  corr_vars,
  use = "pairwise.complete.obs",
  method = "spearman"
)

corr_labels <- c(
  u5mr = "U5MR",
  dtp3 = "DTP3",
  mcv1 = "MCV1",
  pol3 = "Pol3",
  gdp_pc = "GDP pc",
  health_exp_pc = "Health exp. pc",
  fertility = "Fertility",
  basic_water = "Basic water",
  basic_sanitation = "Basic sanitation"
)
dimnames(corr_matrix) <- lapply(dimnames(corr_matrix), function(x) corr_labels[x])

png(
  filename = file.path(figures_dir, "figure6_spearman_correlation_heatmap.png"),
  width = 3000,
  height = 2400,
  res = 400
)

corrplot(
  corr_matrix,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  number.cex = 0.7,
  tl.cex = 0.85,
  tl.col = "black",
  diag = FALSE,
  title = "Spearman correlation matrix",
  mar = c(0, 0, 2, 0)
)

dev.off()

# ============================================================
# Figure 7: GDP vs U5MR scatterplot
# ============================================================

fig7 <- ggplot(immunization_panel, aes(x = gdp_pc, y = u5mr)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "loess", se = TRUE) +
  scale_x_log10(labels = comma) +
  labs(
    title = "Association between GDP per capita and under-five mortality",
    subtitle = "GDP per capita shown on log scale",
    x = "GDP per capita, current US$ (log scale)",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure7_gdp_u5mr_scatter.png"),
  plot = fig7,
  width = 8,
  height = 5,
  dpi = 600
)

# ============================================================
# Figure 8: Fertility vs U5MR scatterplot
# ============================================================

fig8 <- ggplot(immunization_panel, aes(x = fertility, y = u5mr)) +
  geom_point(alpha = 0.45) +
  geom_smooth(method = "loess", se = TRUE) +
  labs(
    title = "Association between fertility rate and under-five mortality",
    subtitle = "Each point represents one country-year observation",
    x = "Total fertility rate, births per woman",
    y = "Under-five mortality rate per 1,000 live births"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(figures_dir, "figure8_fertility_u5mr_scatter.png"),
  plot = fig8,
  width = 8,
  height = 5,
  dpi = 600
)

message("Figures generated successfully.")
