## ======================================
## 0 — working directory  (edit if needed)
## ======================================
#setwd("/Users/hongruzhao/myfile/Research simulation/Federated learning/data_final")
setwd("~/myfile/Research simulation/Federated learning/data_final")

## ======================================
## 1 — helper: load one hospital & append `site`
## ======================================
load_one_site <- function(site, all_sites) {
  X_file <- sprintf("heart_%s_features.csv", site)
  y_file <- sprintf("heart_%s_targets.csv",  site)
  
  if (!file.exists(X_file) || !file.exists(y_file))
    stop("Files not found for site ‘", site, "’")
  
  X <- read.csv(X_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  y <- read.csv(y_file, stringsAsFactors = FALSE)$target
  
  ## add 4-level site factor
  X$site <- factor(site, levels = all_sites)
  
  ## recode categorical predictors (incl. site)
  cat_vars <- c("sex","cp","fbs","restecg","exang","slope","ca","thal","site")
  cat_vars <- intersect(cat_vars, names(X))
  X[cat_vars] <- lapply(X[cat_vars], factor)
  
  list(X = X, y = y)
}

## ======================================
## 2 — load all four hospitals
## ======================================
sites <- c("cleveland", "hungarian", "swiss", "longbeach")
datasets <- lapply(setNames(sites, sites),
                   load_one_site,
                   all_sites = sites)

## grab individual objects if you want them
X_cleveland <- datasets$cleveland$X; y_cleveland <- datasets$cleveland$y
X_hungarian <- datasets$hungarian$X; y_hungarian <- datasets$hungarian$y
X_swiss     <- datasets$swiss$X;     y_swiss     <- datasets$swiss$y
X_longbeach <- datasets$longbeach$X; y_longbeach <- datasets$longbeach$y

## quick check
sapply(datasets, function(d) dim(d$X))   # rows × 14 (13 + site)

## ======================================
## 3 — stack all X, y
## ======================================
X_all <- do.call(rbind, lapply(datasets, `[[`, "X"))
rownames(X_all) <- NULL
y_all <- unlist(lapply(datasets, `[[`, "y"), use.names = FALSE)

## ======================================
## 4 — split variable roles
## ======================================
## Only these are allowed as "global" inputs for testing:
ALLOWED_GLOBALS <- c("age","sex","chol","trestbps","thalach","oldpeak")

## For constructing the *global design block*, we keep cp as in your original
## (so it’s available as a covariate), but it is NOT allowed as a test input.
global_base_vars <- c(ALLOWED_GLOBALS, "cp")

site_vars_lack_missing  <- c("fbs","restecg","exang")
site_vars_huge_missing  <- c("slope","ca","thal")
site_vars               <- c(site_vars_lack_missing, site_vars_huge_missing)

## keep global predictors ***plus site***
X_global <- X_all[ , c(intersect(global_base_vars, names(X_all)), "site") ]

X_site_lack_missing <- X_all[ , site_vars_lack_missing ]
X_site_huge_missing <- X_all[ , site_vars_huge_missing ]

## handle the huge-missing block: cle vs non-cle
X_site_huge_missing_cleveland      <- X_site_huge_missing[X_all$site=="cleveland", ]
X_site_huge_missing_non_cleveland  <- X_site_huge_missing[!(X_all$site=="cleveland"), c("slope","thal")]

## convert to factors and replace NA → "-1"
X_site_huge_missing_non_cleveland[] <- lapply(X_site_huge_missing_non_cleveland, factor)
for (v in names(X_site_huge_missing_non_cleveland)) {
  x <- X_site_huge_missing_non_cleveland[[v]]
  if (!"-1" %in% levels(x)) levels(x) <- c(levels(x), "-1")
  x[is.na(x)] <- "-1"
  X_site_huge_missing_non_cleveland[[v]] <- x
}

## ======================================
## 5 — separate Cleveland vs non-Cleveland
## ======================================
is_cle <- X_all$site == "cleveland"

X_full_cleveland <- cbind(X_global[ is_cle,  ],
                          X_site_lack_missing[ is_cle,  ],
                          X_site_huge_missing_cleveland)

X_full_non_cleveland <- cbind(X_global[!is_cle,  ],
                              X_site_lack_missing[!is_cle, ],
                              X_site_huge_missing_non_cleveland)

## keep flags (complete cases)
keep_cle    <- rowSums(is.na(X_full_cleveland))     == 0
keep_noncle <- rowSums(is.na(X_full_non_cleveland)) == 0

## responses
y_cle    <- y_all[ is_cle ]
y_noncle <- y_all[!is_cle ]

## final design-time data.frames with no NA
X_cle_noNA     <- X_full_cleveland[ keep_cle,  , drop = FALSE ]
y_cle_noNA     <- y_cle[ keep_cle ]

X_noncle_noNA  <- X_full_non_cleveland[ keep_noncle, , drop = FALSE ]
y_noncle_noNA  <- y_noncle[ keep_noncle ]

## ======================================
#################################################################
##  6 — build numeric design matrices with uniform column order
#################################################################

## --------------------------------------------------
## 6.1 GLOBAL block (same for every row)
##     (Build programmatically; no hard-coded indices)
## --------------------------------------------------
global_cols <- intersect(colnames(X_cle_noNA), c("site", global_base_vars))
## Create a combined data frame in the SAME row order as later blocks:
global_df_all <- rbind(X_cle_noNA[, global_cols, drop = FALSE],
                       X_noncle_noNA[, global_cols, drop = FALSE])

form_global <- as.formula(paste("~ ", paste(global_cols, collapse = " + ")))
Z_global_all <- model.matrix(form_global, data = global_df_all)[,-1]

## --------------------------------------------------
## 6.2 SITE-SPECIFIC block (site main-effects + interactions)
## --------------------------------------------------
form_site_cle  <- ~ -1 + (fbs + restecg + exang + slope + ca + thal)
Z_site_cle     <- model.matrix(form_site_cle, data = X_cle_noNA)

form_site_rest <- ~ -1 + site:(fbs + restecg + exang + slope + thal)
Z_site_rest    <- model.matrix(form_site_rest, data = X_noncle_noNA)

## if a column name contains “cleveland”, remove it
clev_mask <- grepl("cleveland", colnames(Z_site_rest), fixed = TRUE)
if (any(clev_mask)) {
  Z_site_rest <- Z_site_rest[ , !clev_mask, drop = FALSE ]
}

## --------------------------------------------------
## 6.3 Align row groups and stack
## --------------------------------------------------
n_cle  <- nrow(Z_site_cle)         
n_rest <- nrow(Z_site_rest)

Z_global_cle  <- Z_global_all[  1:n_cle,  , drop = FALSE ]
Z_global_rest <- Z_global_all[(n_cle+1):nrow(Z_global_all), , drop = FALSE ]

## zero-filled placeholders so col-counts match
zero_cle_block  <- matrix(0, nrow = n_cle,  ncol = ncol(Z_site_rest),
                          dimnames = list(NULL, colnames(Z_site_rest)))
zero_rest_block <- matrix(0, nrow = n_rest, ncol = ncol(Z_site_cle),
                          dimnames = list(NULL, colnames(Z_site_cle)))

## bind columns for each row group
Z_cle_final  <- cbind(Z_global_cle,  Z_site_cle,  zero_cle_block)
Z_rest_final <- cbind(Z_global_rest, zero_rest_block, Z_site_rest)

## one grand design matrix + response
Z_all_final <- rbind(Z_cle_final, Z_rest_final)
y_final     <- c(y_cle_noNA, y_noncle_noNA)

## ======================================
## 7 — Hypothesis test helper
##     Only allows the six pre-approved globals
## ======================================
library(glmtlp)



hetero_LRT <- function(globals,
                       Z = Z_all_final,
                       y = y_final,
                       allowed = ALLOWED_GLOBALS,
                       seed = 2025,
                       family = "binomial",
                       penalty = "l0",
                       nfolds = 5) {
  ## 1) Validate and normalize input
  if (missing(globals) || length(globals) == 0)
    stop("Please supply at least one variable to test.")
  globals <- unique(tolower(as.character(globals)))
  
  bad <- setdiff(globals, tolower(allowed))
  # if (length(bad) > 0) {
  #   stop("Invalid global(s): ", paste(bad, collapse = ", "),
  #        ". Allowed: ", paste(allowed, collapse = ", "), call. = FALSE)
  # }
  
  ## 2) Identify GLOBAL block columns to test (prefix match)
  pat <- paste0("^(", paste(globals, collapse = "|"), ")")
  B_index <- grep(pat, colnames(Z))
  if (length(B_index) == 0)
    stop("No matching columns in the design matrix for: ",
         paste(globals, collapse = ", "), call. = FALSE)
  
  ## 3) Build H0 / H1 design matrices
  Z_H1 <- Z
  Z_H0 <- Z[, -B_index, drop = FALSE]
  
  ## 4) Cross-validated fits (H0 & H1)
  set.seed(seed)
  
  ## H0: standard cv fit
  cv_H0 <- cv.glmtlp(X = Z_H0, y = y, family = family,
                     penalty = penalty, nfolds = nfolds)
  
  ## H1: unpenalize the global block by setting its penalty.factor to 0
  p  <- ncol(Z_H1)
  pf <- rep(1, p); pf[B_index] <- 0
  cv_H1 <- cv.glmtlp(X = Z_H1, y = y, family = family,
                     penalty = penalty, penalty.factor = pf, nfolds = nfolds)
  
  ## 5) Extract coefficients at λ_min (keep NAMES!)
  c0 <- coef(cv_H0)  # named vector: "(Intercept)", terms...
  c1 <- coef(cv_H1)
  
  ## 6) Compute fitted linear predictors (η = Xβ + intercept)
  b0 <- as.numeric(c0)
  b1 <- as.numeric(c1)
  eta0 <- as.vector(Z_H0 %*% b0[-1] + b0[1])
  eta1 <- as.vector(Z_H1 %*% b1[-1] + b1[1])
  
  ## 7) Binomial log-likelihoods
  ll0 <- sum(y * eta0 - log1p(exp(eta0)))
  ll1 <- sum(y * eta1 - log1p(exp(eta1)))
  
  ## 8) LRT
  LR_stat <- -2 * (ll0 - ll1)
  df      <- length(B_index)
  p_val   <- pchisq(LR_stat, df = df, lower.tail = FALSE)
  
  ## NEW: supports at λ_min (exactly as requested)
  nzcoef_H0   <- c0[c0 != 0]
  nzcoef_H1   <- c1[c1 != 0]
  support_H0  <- setdiff(names(nzcoef_H0), "(Intercept)")
  support_H1  <- setdiff(names(nzcoef_H1), "(Intercept)")
  added_in_H1   <- setdiff(support_H1, support_H0)
  dropped_in_H1 <- setdiff(support_H0, support_H1)
  
  out <- list(
    tested_globals = globals,
    df             = df,
    LR             = LR_stat,
    p_value        = p_val,
    B_index        = B_index,
    B_columns      = colnames(Z)[B_index],
    cv_H0          = cv_H0,
    cv_H1          = cv_H1,
    ## NEW: supports & nonzero coefs
    nzcoef_H0      = nzcoef_H0,
    nzcoef_H1      = nzcoef_H1,
    support_H0     = support_H0,
    support_H1     = support_H1,
    added_in_H1    = added_in_H1,
    dropped_in_H1  = dropped_in_H1
  )
  class(out) <- c("hetero_LRT_result", class(out))
  out
}


print.hetero_LRT_result <- function(x, max_terms = 30, ...) {
  cat("Heterogeneous LR test\n")
  cat("  Tested globals :", paste(x$tested_globals, collapse = ", "), "\n")
  cat("  df             :", x$df, "\n")
  cat("  LR statistic   :", signif(x$LR, 5), "\n")
  cat("  p-value        :", signif(x$p_value, 5), "\n")
  
  if (!is.null(x$B_columns)) {
    cat("  Columns tested :", paste(x$B_columns, collapse = ", "), "\n")
  }
  
  # helper for pretty truncation
  prettify <- function(v, k = max_terms) {
    if (length(v) == 0) return("<none>")
    if (length(v) <= k)  return(paste(v, collapse = ", "))
    paste(paste(v[seq_len(k)], collapse = ", "),
          sprintf("... (+%d more)", length(v) - k))
  }
  
  
  # Optionally, show nonzero coefficient vectors (incl. intercept) in compact form
  cat("\n  Nonzero coefficients under H0 (incl. intercept):\n")
  print(x$nzcoef_H0)
  cat("\n  Nonzero coefficients under H1 (incl. intercept):\n")
  print(x$nzcoef_H1)
  
  invisible(x)
}

## ======================================
## 8 — Examples (uncomment to run)
## ======================================
res_age  <- hetero_LRT("age")
res_sex  <- hetero_LRT("sex")
res_chol <- hetero_LRT("chol")
res_cp <- hetero_LRT("cp")

print(res_age)
print(res_sex)
print(res_chol)
print(res_cp)

# ## test multiple globals jointly
# res_two <- hetero_LRT(c("age","sex"))
# print(res_two)

# ## run all six and collect p-values


results_all <- lapply(ALLOWED_GLOBALS, function(v) {
   r <- hetero_LRT(v)
   data.frame(var = v, df = r$df, LR = r$LR, p = r$p_value)
 })
do.call(rbind, results_all)





ALLOWED_TOKENS <- c(
  # global (as before)
  "age","sex","chol","trestbps","thalach","oldpeak",
  # optional global you had kept as covariate
  "cp",
  # site main-effect dummies in Z_global_all (e.g., sitehungarian, siteswiss, ...)
  "site"#,
  # site-specific blocks
#  "fbs","restecg","exang","slope","ca","thal"#,
  # non-Cleveland site-specific interactions only:
#  "site:fbs","site:restecg","site:exang","site:slope","site:thal"
)







results_all <- lapply(ALLOWED_TOKENS, function(v) {
  r <- hetero_LRT(v)
  data.frame(var = v, df = r$df, LR = r$LR, p = r$p_value)
})
do.call(rbind, results_all)

