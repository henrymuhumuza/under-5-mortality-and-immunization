# Response to Reviewer and Editor

**Submission ID:** b2bd4114-5629-44f6-918a-bd9e61c4abe1  
**Manuscript:** *Association Between Routine Childhood Immunization Coverage and Under-Five Mortality in Sub-Saharan Africa: An Ecological Panel Study, 2000–2024*

We thank the reviewer and editor for the constructive recommendations. We have revised the manuscript and reproducible analysis workflow. Throughout, results are framed as country-level associations, not causal effects.

| Comment | Response |
|---|---|
| Account for serial correlation and cross-sectional dependence. | We retained country-clustered standard errors and added Driscoll–Kraay standard errors, which allow serial and cross-sectional dependence, plus two-way clustered standard errors. These results appear in new Supplementary Table S6. |
| Justify linear modeling and add robustness checks. | We clarify that U5MR is continuous and the linear scale gives directly interpretable differences in deaths per 1,000 live births. We added a quadratic DTP3 specification (Supplementary Table S8), country-specific linear trends, and dependence-robust inference. The quadratic diagnostic is reported as a robustness check, not as causal evidence. |
| Examine alternative lags and country-specific trends. | We added prospective 1-, 2-, and 3-year lag models with Driscoll–Kraay inference (Supplementary Table S7), and a country-specific linear-trend model (Supplementary Table S6). |
| Address high antigen correlations and VIFs. | Antigens remain estimated in separate models. The Methods, Results, and Discussion now explicitly state that high correlations and VIFs make simultaneous antigen coefficients unstable and unsuitable for antigen-specific causal comparisons. |
| Examine complete-case analysis and missing 2024 health-expenditure data. | We added an included-versus-excluded comparison (Supplementary Table S5) and matched analyses omitting health expenditure in 2000–2023 and 2000–2024 (Supplementary Table S6). The manuscript now distinguishes the assembled 2000–2024 panel from the 2000–2023 primary complete-case sample. |
| Avoid causal language. | We replaced causal terminology with “associated with” wording and expanded limitations on residual confounding, reverse causation, measurement error, and ecological fallacy. |
| Edit tables, figures, and references. | Table/figure notes were clarified and new Supplementary Tables S5–S8 were added. All additions are reproducible in `scripts/10_revision_robustness.R`. |

We believe these revisions strengthen the methodological transparency and the appropriately cautious interpretation of the study.
