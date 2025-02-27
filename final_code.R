#############################################
# tlp.R: Global Model using glmtlp          
# (No cohort-specific effects or response standardization)
#############################################

# -------------------------------
# 1. Data Loading and Cleaning
# -------------------------------
# Set the working directory to where the dataset is stored.
setwd("C:/Research in R/Cancer")

# Read the METABRIC RNA Mutation CSV file into a data frame.
# header = TRUE indicates that the first row contains column names.
# stringsAsFactors = FALSE prevents character strings from being automatically converted to factors.
metabric_mutation <- read.csv("METABRIC_RNA_Mutation.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove columns that contain any NA values.
# sapply() applies a function over each column to check for NAs.
# We keep only columns where no NA is found.
metabric_mutation_clean <- metabric_mutation[, !sapply(metabric_mutation, function(x) any(is.na(x)))]

# Print the dimensions of the data before and after cleaning for verification.
cat("Dimensions before cleaning:", dim(metabric_mutation), "\n")
cat("Dimensions after cleaning:", dim(metabric_mutation_clean), "\n")

# Remove _mut columns with more than d unique levels.
# Set threshold for maximum allowed unique levels.
d <- 400
# Identify column names that end with '_mut'
mut_cols <- grep("_mut$", names(metabric_mutation_clean), value = TRUE)
# Count the number of unique values for each _mut column.
mut_levels <- sapply(metabric_mutation_clean[, mut_cols], function(x) length(unique(x)))
# Determine which _mut columns exceed the threshold.
mut_to_remove <- names(mut_levels)[mut_levels > d]
cat("Removing _mut columns with more than", d, "levels:\n")
print(mut_to_remove)
# Remove those columns from the dataset.
metabric_mutation_clean2 <- metabric_mutation_clean[, !(names(metabric_mutation_clean) %in% mut_to_remove)]

# -------------------------------
# 2. Define Predictor Groups and Build Design Matrices
# -------------------------------
# Specify columns that should not be used as predictors.
non_predictors <- c("patient_id", "overall_survival_months", "overall_survival", "death_from_cancer")

# Define global clinical predictors.
global_clinical <- c("age_at_diagnosis", "cancer_type", "cancer_type_detailed", 
                     "pam50_._claudin.low_subtype", "er_status_measured_by_ihc", "er_status", 
                     "her2_status_measured_by_snp6", "her2_status", "inferred_menopausal_state", 
                     "integrative_cluster", "primary_tumor_laterality", "lymph_nodes_examined_positive", 
                     "nottingham_prognostic_index", "oncotree_code", "pr_status", 
                     "X3.gene_classifier_subtype", "cohort")

# Define site-specific clinical predictors.
site_specific_clinical <- c("type_of_breast_surgery", "cellularity", "chemotherapy", 
                            "tumor_other_histologic_subtype", "hormone_therapy", "radio_therapy")

# Combine all columns that should be excluded from gene predictors.
all_excluded <- c(non_predictors, global_clinical, site_specific_clinical)
# Identify gene predictors by removing all excluded columns from the dataset.
gene_predictors <- setdiff(names(metabric_mutation_clean2), all_excluded)

# Define global predictors as a combination of global clinical variables and gene predictors.
global_predictors <- c(global_clinical, gene_predictors)
# Site-specific predictors remain as defined.
site_specific_predictors <- site_specific_clinical

# Print the predictor groups to verify.
cat("\nGlobal Predictors:\n")
print(global_predictors)
cat("\nSite-Specific Predictors:\n")
print(site_specific_predictors)

# Build the design matrices for modeling.
# (A) Global design matrix: subset the data to only include global predictors.
global_data <- metabric_mutation_clean2[, global_predictors]
# Use model.matrix() to create a design matrix; "~ . - 1" includes all predictors and omits the intercept.
global_design <- model.matrix(~ . - 1, data = global_data)

# (B) Site-specific design matrix: create interaction terms between 'cohort' and each site-specific predictor.
# The formula includes "0 +" to omit the intercept.
interaction_formula <- as.formula(paste("~ 0 + cohort:(", paste(site_specific_predictors, collapse = " + "), ")"))
site_specific_design <- model.matrix(interaction_formula, data = metabric_mutation_clean2)

# (C) Full design matrix: combine the global and site-specific design matrices.
full_design <- cbind(global_design, site_specific_design)

# -------------------------------
# 3. Define the Response
# -------------------------------
# Use the log-transformed overall survival months as the response variable.
# Adding 1 inside log() helps avoid issues if any survival value is zero.
y_log <- log(metabric_mutation_clean2$overall_survival_months + 1)

# -------------------------------
# 4. Global Model: 1000 Iterations for Train-Test Split and Model Evaluation
# -------------------------------
# Load the glmtlp package for fitting penalized regression models.
library(glmtlp)

# Set the number of iterations for repeated train-test splits.
nIter <- 1000
# Initialize a vector to store test Mean Squared Error (MSE) for each iteration.
test_mse_all <- numeric(nIter)
set.seed(2025)  # Ensure reproducibility of random splits

# Determine the total number of observations.
n <- nrow(full_design)
for (iter in 1:nIter) {
  # Randomly split the data: 70% for training and 30% for testing.
  train_idx <- sample(seq_len(n), size = floor(0.7 * n))
  test_idx <- setdiff(seq_len(n), train_idx)
  
  # Create training and test sets for predictors and response.
  X_train <- full_design[train_idx, , drop = FALSE]
  y_train <- y_log[train_idx]
  
  X_test <- full_design[test_idx, , drop = FALSE]
  y_test <- y_log[test_idx]
  
  # Set penalty factors for each predictor (here all predictors are penalized equally).
  penalty_factor <- rep(1, ncol(X_train))
  
  # Fit the model using 5-fold cross-validation with the tlp penalty.
  cv_fit <- cv.glmtlp(
    X_train, y_train,
    type = "ls",       # Using least squares (Gaussian) loss function.
    penalty = "tlp",
    penalty.factor = penalty_factor,
    nfolds = 5
  )
  
  # Predict the response on the test set.
  y_test_pred <- predict(cv_fit, X_test)
  # Calculate residuals and compute the mean squared error for the test set.
  resid_test <- y_test - y_test_pred
  test_mse_all[iter] <- mean(resid_test^2)
  
  # Print progress every 100 iterations along with the current average MSE.
  if(iter %% 100 == 0) {
    cat("Iteration", iter, "completed. Current average MSE:", mean(test_mse_all[1:iter]), "\n")
  }
}

# Compute and print the overall average test MSE across all iterations.
avg_test_mse <- mean(test_mse_all)
cat("\nAverage Test MSE (on log-scale) over", nIter, "iterations =", avg_test_mse, "\n")

# Write the test MSE values for all iterations to a CSV file.
write.csv(test_mse_all, file = "test_mse_all.csv", row.names = FALSE)

#############################################
# sep.R: Separate CV.glmtlp by Cohort       
# for METABRIC_RNA_Mutation dataset         
#############################################

# ----- Data Loading and Cleaning -----
# Set the working directory and load the dataset.
setwd("C:/Research in R/Cancer")
metabric_mutation <- read.csv("METABRIC_RNA_Mutation.csv", header = TRUE, stringsAsFactors = FALSE)

# Remove columns that contain any NA values.
metabric_mutation_clean <- metabric_mutation[, !sapply(metabric_mutation, function(x) any(is.na(x)))]
cat("Dimensions before cleaning:", dim(metabric_mutation), "\n")
cat("Dimensions after cleaning:", dim(metabric_mutation_clean), "\n")

# Remove _mut columns with more than a specified number of unique levels.
d <- 400
mut_cols <- grep("_mut$", names(metabric_mutation_clean), value = TRUE)
mut_levels <- sapply(metabric_mutation_clean[, mut_cols], function(x) length(unique(x)))
mut_to_remove <- names(mut_levels)[mut_levels > d]
cat("Removing _mut columns with more than", d, "levels:\n")
print(mut_to_remove)
metabric_mutation_clean2 <- metabric_mutation_clean[, !(names(metabric_mutation_clean) %in% mut_to_remove)]

# ----- Define Predictor Groups -----
# Specify non-predictor columns that should be excluded.
non_predictors <- c("patient_id", "overall_survival_months", "overall_survival", "death_from_cancer")

# Define global clinical predictors.
global_clinical <- c("age_at_diagnosis", "cancer_type", "cancer_type_detailed", 
                     "pam50_._claudin.low_subtype", "er_status_measured_by_ihc", "er_status", 
                     "her2_status_measured_by_snp6", "her2_status", "inferred_menopausal_state", 
                     "integrative_cluster", "primary_tumor_laterality", "lymph_nodes_examined_positive", 
                     "nottingham_prognostic_index", "oncotree_code", "pr_status", 
                     "X3.gene_classifier_subtype", "cohort")

# Define site-specific clinical predictors.
site_specific_clinical <- c("type_of_breast_surgery", "cellularity", "chemotherapy", 
                            "tumor_other_histologic_subtype", "hormone_therapy", "radio_therapy")

# Exclude non-predictors, global clinical, and site-specific clinical predictors to isolate gene predictors.
all_excluded <- c(non_predictors, global_clinical, site_specific_clinical)
gene_predictors <- setdiff(names(metabric_mutation_clean2), all_excluded)

# Combine global clinical and gene predictors to form the global predictor set.
global_predictors <- c(global_clinical, gene_predictors)
# Site-specific predictors remain as defined.
site_specific_predictors <- site_specific_clinical

# Print the predictor groups to verify.
cat("\nGlobal Predictors:\n")
print(global_predictors)
cat("\nSite-Specific Predictors:\n")
print(site_specific_predictors)

# ----- Build the Design Matrices -----
# Global design matrix: subset the data to the global predictor columns and create a model matrix (without intercept).
global_data <- metabric_mutation_clean2[, global_predictors]
global_design <- model.matrix(~ . - 1, data = global_data)

# Site-specific design matrix: build a model matrix for site-specific predictors.
# Here we simply include the site-specific predictors (without interacting with cohort).
# (If interaction with cohort was desired, you could modify the formula accordingly.)
interaction_formula <- as.formula(paste("~ 0 +(", paste(site_specific_predictors, collapse = " + "), ")"))
site_specific_design <- model.matrix(interaction_formula, data = metabric_mutation_clean2)

# Full design matrix: combine the global and site-specific design matrices.
full_design <- cbind(global_design, site_specific_design)

# ----- Define the Response Variables -----
# Define the response variable as the log-transformed overall survival months.
# Adding 1 ensures that we do not take the log of zero.
y_log <- log(metabric_mutation_clean2$overall_survival_months + 1)

##############################################
#  Separate CV.glmtlp for Each Cohort (1000 iterations)  
##############################################
# Load the glmtlp package for model fitting.
library(glmtlp)

# Identify the unique cohorts in the dataset.
cohorts <- unique(metabric_mutation_clean2$cohort)
nCohorts <- length(cohorts)

# Initialize a matrix to store the test MSE for each cohort across iterations.
nIter <- 1000
mse_mat <- matrix(NA, nrow = nIter, ncol = nCohorts)
colnames(mse_mat) <- cohorts

set.seed(2025)  # For reproducibility

# Loop over the specified number of iterations.
for (iter in 1:nIter) {
  cat("Iteration", iter, "\n")
  
  # Loop over each cohort.
  for (i in seq_along(cohorts)) {
    cohort <- cohorts[i]
    # Subset the data for the current cohort.
    idx_cohort <- metabric_mutation_clean2$cohort == cohort
    X_cohort <- full_design[idx_cohort, , drop = FALSE]
    y_cohort <- y_log[idx_cohort]
    
    # Determine the number of observations for the current cohort.
    n_cohort <- nrow(X_cohort)
    if(n_cohort < 5) next  # Skip if there are too few observations to split.
    
    # Perform a random 70/30 train-test split within the current cohort.
    train_idx <- sample(seq_len(n_cohort), size = floor(0.7 * n_cohort))
    test_idx <- setdiff(seq_len(n_cohort), train_idx)
    
    # Create training and test sets.
    X_train <- X_cohort[train_idx, , drop = FALSE]
    y_train <- y_cohort[train_idx]
    X_test <- X_cohort[test_idx, , drop = FALSE]
    y_test <- y_cohort[test_idx]
    
    # Set penalty factors so that all predictors are equally penalized.
    penalty_factor <- rep(1, ncol(X_train))
    
    # Fit the model using 5-fold cross-validation with the tlp penalty.
    cv_fit <- cv.glmtlp(
      X_train, y_train,
      type = "ls",         # Least squares loss (Gaussian)
      penalty = "tlp",
      penalty.factor = penalty_factor,
      nfolds = 5
    )
    
    # Predict responses on the test set.
    pred_test <- predict(cv_fit, X_test)
    
    # Calculate the test MSE for the current cohort and iteration.
    mse_test <- mean((y_test - pred_test)^2)
    mse_mat[iter, cohort] <- mse_test
  }
}

# Compute the average test MSE for each cohort over all iterations.
avg_test_mse <- colMeans(mse_mat, na.rm = TRUE)
cat("\nAverage Test MSE by Cohort over", nIter, "iterations:\n")
print(avg_test_mse)

# ----- Determine the Number of Test Observations per Cohort -----
# Initialize a vector to store the number of test observations for each cohort.
n_valid_vec <- numeric(length(cohorts))
names(n_valid_vec) <- cohorts

# For each cohort, perform a train-test split and record the number of test observations.
for (i in seq_along(cohorts)) {
  cohort <- cohorts[i]
  idx_cohort <- metabric_mutation_clean2$cohort == cohort
  X_cohort <- full_design[idx_cohort, , drop = FALSE]
  y_cohort <- y_log[idx_cohort]
  
  n_cohort <- nrow(X_cohort)
  if(n_cohort < 5) {
    n_valid_vec[i] <- NA  # Mark as NA if too few observations.
    next
  }
  
  train_idx <- sample(seq_len(n_cohort), size = floor(0.7 * n_cohort))
  test_idx <- setdiff(seq_len(n_cohort), train_idx)
  n_valid_vec[i] <- length(test_idx)
}

# Print the number of test observations per cohort.
print(n_valid_vec)

# ----- Compute Overall Weighted Test MSE -----
# Create a weight matrix by replicating the number of test observations for each cohort across iterations.
weight_matrix <- matrix(n_valid_vec, nrow = nIter, ncol = length(n_valid_vec), byrow = TRUE)

# For each iteration, compute the weighted sum of the MSE values.
weighted_sum_per_iter <- rowSums(mse_mat * weight_matrix)

# Compute the weighted average MSE for each iteration.
weighted_avg_per_iter <- weighted_sum_per_iter / sum(n_valid_vec)

# Compute the overall weighted average test MSE across all iterations.
overall_weighted_avg <- mean(weighted_avg_per_iter)
cat("Overall Weighted Average Test MSE =", overall_weighted_avg, "\n")

# Print the mean and standard deviation of the weighted average MSE across iterations.
mean(weighted_avg_per_iter)
sd(weighted_avg_per_iter)

# Read in the overall test MSE from the global model (previously saved in "test_mse_all.csv").
test_mse_df_full <- read.csv("test_mse_all.csv", header = TRUE, stringsAsFactors = FALSE)
mean(test_mse_df_full$x)
sd(test_mse_df_full$x)

# Create a named vector that summarizes the four key metrics:
# mean and standard deviation of the weighted MSE (separated by cohort) 
# and mean and standard deviation of the test MSE from the global model.
result_vector <- c(
  mean_weighted = mean(weighted_avg_per_iter),
  sd_weighted   = sd(weighted_avg_per_iter),
  mean_test     = mean(test_mse_df_full$x),
  sd_test       = sd(test_mse_df_full$x)
)

# Convert the result vector into a one-row data frame with descriptive column names.
result_df <- data.frame(t(result_vector))
colnames(result_df) <- c("mean_sep", "sd_sep", "mean_full", "sd_full")
print(result_df)

