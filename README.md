 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a/README.md b/README.md
new file mode 100644
index 0000000000000000000000000000000000000000..aa429573e698b9fc76629ce6fcdd103f45ef72b7
--- /dev/null
+++ b/README.md
@@ -0,0 +1,54 @@
+# Distributed Algorithms for High-Dimensional Statistical Inference and Structure Learning with Heterogeneous Data
+
+**Hongru Zhao** and **Xiaotong Shen**  
+*School of Statistics, University of Minnesota, Twin Cities*
+
+## Overview
+This repository accompanies the research paper *Distributed Algorithms for High-Dimensional Statistical Inference and Structure Learning with Heterogeneous Data*. It gathers simulation code, heterogeneous inference routines, and the hospital-specific heart disease data used to benchmark the proposed distributed likelihood-ratio tests. The R scripts focus on sparse model selection with \(\ell_0\) penalties, cross-validation, and post-selection inference in both Gaussian and logistic settings.
+
+## Repository structure
+| Path | Description |
+| --- | --- |
+| `ROC_LM.R` | Monte Carlo pipeline for Gaussian responses that generates correlated multi-site designs, fits \(\ell_0\)-penalized linear models, and summarizes full ROC curves alongside the validation-selected point.【F:ROC_LM.R†L1-L171】 |
+| `ROC_logistic.R` | Logistic analogue of the ROC study that tunes \(\ell_0\)-penalized classifiers by cross-entropy and records ROC trajectories for distributed sites.【F:ROC_logistic.R†L1-L203】 |
+| `power_final_high_rho_paper_version.R` | High-dimensional power study that compares heterogeneous likelihood-ratio tests against lasso-based benchmarks over grids of signal strengths, block sizes, and dimensionalities.【F:power_final_high_rho_paper_version.R†L1-L224】 |
+| `heart_disease_two_step_test_paper_version.R` | Two-step heterogeneous logistic testing workflow: builds cross-entropy CV for `glmtlp`, engineers global and site-specific feature blocks from four hospitals, and runs the combined LR test with detailed reporting utilities.【F:heart_disease_two_step_test_paper_version.R†L1-L523】 |
+| `heart_*.csv` | Hospital-level covariate matrices (`*_features.csv`) and binary outcomes (`*_targets.csv`) for Cleveland, Hungarian, Swiss, and Long Beach cohorts used in the heterogeneous testing case study.【F:heart_disease_two_step_test_paper_version.R†L199-L291】 |
+
+## Requirements
+All scripts are written for R (\>= 4.0). Install the following packages before running the code:
+
+```r
+install.packages(c("MASS", "glmnet"))
+# glmtlp is available from GitHub: remotes::install_github("ChongWu-Biostat/glmtlp")
+```
+
+Each script loads only the packages it needs—`MASS` and `glmtlp` for the ROC experiments, and `MASS`, `glmtlp`, plus `glmnet` for the power analyses.【F:ROC_LM.R†L3-L5】【F:ROC_logistic.R†L9-L11】【F:power_final_high_rho_paper_version.R†L1-L7】 The heterogeneous testing workflow also depends solely on `glmtlp` after data preparation.【F:heart_disease_two_step_test_paper_version.R†L1-L523】
+
+## Reproducing the simulations
+1. **Linear-model ROC curves.** Run `ROC_LM.R` in R to simulate multi-site Gaussian data, fit \(\ell_0\)-penalized regressions, and save `roc_plots_LM.pdf` with the aggregated ROC curves.【F:ROC_LM.R†L62-L171】 Adjust the `params1`/`params2` lists to explore additional designs.【F:ROC_LM.R†L127-L146】
+2. **Logistic ROC curves.** Execute `ROC_logistic.R` to generate Bernoulli responses under heterogeneous designs, tune by validation cross-entropy, and write `roc_logistic_ce.pdf` for both parameter regimes.【F:ROC_logistic.R†L75-L203】
+3. **Power study under high correlation.** Launch `power_final_high_rho_paper_version.R` to reproduce the Monte Carlo experiment that contrasts the proposed heterogeneous LR test against a lasso-split benchmark. The script iterates over grids of dimensionalities, block sizes, and signal strengths, collecting rejection rates into `results_df_all`.【F:power_final_high_rho_paper_version.R†L56-L224】 Large replication counts (`n_rep = 1000`) may require parallelization for faster turnaround.【F:power_final_high_rho_paper_version.R†L63-L105】
+
+## Heterogeneous heart disease analysis
+The heart disease workflow illustrates how to apply the distributed test with real-world heterogeneity:
+
+1. **Data assembly.** `load_one_site()` reads each hospital’s features and outcomes, adds a site factor, and recodes categorical predictors.【F:heart_disease_two_step_test_paper_version.R†L199-L217】 The datasets from Cleveland, Hungarian, Swiss, and Long Beach are concatenated to form a pooled sample.【F:heart_disease_two_step_test_paper_version.R†L223-L236】
+2. **Feature engineering.** Global covariates shared by all sites, site-level main effects, and interaction terms are constructed with consistent ordering before stacking into a block-diagonal design matrix that preserves heterogeneity.【F:heart_disease_two_step_test_paper_version.R†L241-L340】 Missingness in heavily incomplete site-specific variables is imputed with sentinel levels to retain patients from non-Cleveland hospitals.【F:heart_disease_two_step_test_paper_version.R†L252-L263】
+3. **Testing routine.** `hetero_LRT()` cross-validates an \(\ell_0\)-penalized null model, fits a logistic GLM on the union support, and reports likelihood-ratio statistics, p-values, and selected supports with a dedicated print method for inspection.【F:heart_disease_two_step_test_paper_version.R†L356-L523】 Example calls at the end of the script showcase single-variable and joint global tests.【F:heart_disease_two_step_test_paper_version.R†L504-L523】
+
+To reproduce the paper’s heart disease results, adjust the working directory to the folder containing the CSV files (see the commented `setwd()` lines near the top of the script) and run the driver examples or loop over the allowed global variables.【F:heart_disease_two_step_test_paper_version.R†L190-L205】【F:heart_disease_two_step_test_paper_version.R†L504-L523】 The script returns structured objects so you can tabulate or visualize outcomes across hospitals.
+
+## Outputs
+Each simulation script emits publication-ready artifacts:
+
+- `roc_plots_LM.pdf` and `roc_logistic_ce.pdf` provide side-by-side ROC panels for two representative parameter settings.【F:ROC_LM.R†L148-L171】【F:ROC_logistic.R†L169-L203】
+- `results_df_all` from the power study records empirical rejection rates for the heterogeneous LR test, its normal approximation, and the lasso-split benchmark, facilitating downstream plotting or table generation.【F:power_final_high_rho_paper_version.R†L199-L224】
+- `hetero_LRT()` returns a rich summary that the custom `print()` method formats for reports or appendices, including nonzero coefficients under both hypotheses.【F:heart_disease_two_step_test_paper_version.R†L439-L523】
+
+## Citation
+If you use this repository, please cite the accompanying paper:
+
+> Hongru Zhao and Xiaotong Shen. *Distributed Algorithms for High-Dimensional Statistical Inference and Structure Learning with Heterogeneous Data*. School of Statistics, University of Minnesota, Twin Cities.
+
+We appreciate references to the scripts or datasets when reusing them in derivative work. 
EOF
)
