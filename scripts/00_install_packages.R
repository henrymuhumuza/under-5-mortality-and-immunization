required_packages <- c(
  "tidyverse",
  "readxl",
  "janitor",
  "psych",
  "corrplot",
  "GGally",
  "patchwork",
  "countrycode",
  "broom",
  "modelsummary",
  "gtsummary",
  "gt",
  "flextable",
  "officer",
  "fixest",
  "plm",
  "performance",
  "car",
  "forestplot",
  "ggforestplot"
)

new_packages <- required_packages[
  !(required_packages %in% installed.packages()[,"Package"])
]

if(length(new_packages)) {
  install.packages(new_packages)
}

invisible(lapply(required_packages, library, character.only = TRUE))

message("All packages loaded successfully.")