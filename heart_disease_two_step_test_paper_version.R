## ===============================================================
##  Packages
## ===============================================================
library(glmtlp)

## ===============================================================
##  CV function for glmtlp (binomial, L0) using cross-entropy
##  One-dimensional path index over paired (lambda_i, kappa_i)
## ===============================================================
# Returns:
#   $pairs      : data.frame(idx, lambda, kappa)
#   $idx_min    : selected index minimizing total CV CE
#   $fold_loss  : nfolds x L matrix of CE per fold & index
#   $cv_total   : length-L vector (sum over folds)
#   $cv_mean    : cv_total / n
#   $fit        : final full-data glmtlp() object along the same path
#   S3: coef(cv_obj) returns coefficients at the selected index
cv_glmtlp_ce <- function(
    X, y,
    nfolds         = 5,
    foldid         = NULL,   # optional vector in {1,..,nfolds}
    seed           = NULL,
    penalty.factor = NULL,
    standardize    = FALSE,
    lambda         = NULL,   # optional: enforce this path
    kappa          = NULL,   # optional: enforce this path
    ...
){
  if (!is.matrix(X)) X <- as.matrix(X)
  storage.mode(X) <- "double"
  n <- nrow(X); p <- ncol(X)
  if (length(y) != n) stop("X and y have incompatible sizes")
  
  # --- binomial 0/1 response
  if (is.factor(y)) {
    if (nlevels(y) != 2) stop("y must be binary (2-level factor).")
    y <- as.integer(y) - 1L
  }
  y <- as.numeric(y)
  if (!all(y %in% c(0, 1))) stop("y must be in {0,1} for binomial.")
  
  # --- penalty factors
  if (is.null(penalty.factor)) penalty.factor <- rep(1, p)
  if (length(penalty.factor) != p) stop("penalty.factor must have length ncol(X).")
  
  # --- folds (stratified by class, unless provided)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(foldid)) {
    k <- nfolds
    pos <- which(y == 1); neg <- which(y == 0)
    foldid <- integer(n)
    foldid[pos] <- sample(rep(1:k, length.out = length(pos)))
    foldid[neg] <- sample(rep(1:k, length.out = length(neg)))
  } else {
    foldid <- as.integer(foldid)
    if (length(foldid) != n) stop("foldid must have length nrow(X).")
    if (!all(foldid %in% seq_len(nfolds))) stop("foldid values must be in {1,..,nfolds}.")
  }
  
  # --- Master fit to get a common path (lambda/kappa)
  master_fit <- glmtlp(
    X, y,
    family          = "binomial",
    penalty         = "l0",
    standardize     = standardize,
    penalty.factor  = penalty.factor,
    lambda          = lambda,
    kappa           = kappa,
    ...
  )
  
  # Extract path sequences
  lambda_seq <- tryCatch(master_fit$lambda, error = function(e) NULL)
  kappa_seq  <- tryCatch(master_fit$kappa,  error = function(e) NULL)
  
  # Number of solutions we can evaluate (prefer intercept length if present)
  L <- tryCatch(length(master_fit$intercept), error = function(e) NULL)
  if (is.null(L) || L < 1) {
    L <- if (!is.null(lambda_seq)) length(lambda_seq) else length(kappa_seq)
  }
  if (is.null(L) || L < 1) stop("glmtlp() did not return a valid penalty path.")
  
  lambda_for_pairs <- rep(NA_real_, L)
  if (!is.null(lambda_seq)) {
    idx_lambda <- seq_len(min(L, length(lambda_seq)))
    lambda_for_pairs[idx_lambda] <- as.numeric(lambda_seq)[idx_lambda]
  }
  kappa_for_pairs <- rep(NA_real_, L)
  if (!is.null(kappa_seq)) {
    idx_kappa <- seq_len(min(L, length(kappa_seq)))
    kappa_for_pairs[idx_kappa] <- as.numeric(kappa_seq)[idx_kappa]
  }
  
  pairs <- data.frame(
    idx    = seq_len(L),
    lambda = lambda_for_pairs,
    kappa  = kappa_for_pairs
  )
  
  # --- 5-fold cross-entropy over fixed path
  fold_loss <- matrix(NA_real_, nrow = nfolds, ncol = L,
                      dimnames = list(paste0("fold", 1:nfolds), paste0("idx", 1:L)))
  for (k in seq_len(nfolds)) {
    test_id  <- which(foldid == k)
    train_id <- setdiff(seq_len(n), test_id)
    
    fit_k <- glmtlp(
      X[train_id, , drop = FALSE], y[train_id],
      family          = "binomial",
      penalty         = "l0",
      standardize     = standardize,
      penalty.factor  = penalty.factor,
      lambda          = lambda_seq,
      kappa           = kappa_seq,
      ...
    )
    
    # Predict probabilities along the entire path on the test fold
    P <- tryCatch(
      predict(fit_k, X    = X[test_id, , drop = FALSE],
              type = "response", which = seq_len(L)),
      error = function(e) {
        tryCatch(
          predict(fit_k, newx = X[test_id, , drop = FALSE],
                  type = "response", which = seq_len(L)),
          error = function(e2) {
            predict(fit_k, newx = X[test_id, , drop = FALSE],
                    type = "response", s = seq_len(L))
          }
        )
      }
    )
    P <- as.matrix(P)
    if (ncol(P) != L) {
      if (length(P) == length(test_id) * L) dim(P) <- c(length(test_id), L)
      else stop("Could not obtain predictions for the full path.")
    }
    
    # Cross-entropy per index
    eps <- 1e-15
    P <- pmin(pmax(P, eps), 1 - eps)
    yte <- matrix(y[test_id], nrow = length(test_id), ncol = L)
    ce  <- colSums(-(yte * log(P) + (1 - yte) * log(1 - P)))
    fold_loss[k, ] <- ce
  }
  
  cv_total <- colSums(fold_loss)
  cv_mean  <- cv_total / n
  idx_min  <- which.min(cv_total)
  
  # --- Final full-data fit along the same path
  final_fit <- glmtlp(
    X, y,
    family          = "binomial",
    penalty         = "l0",
    standardize     = standardize,
    penalty.factor  = penalty.factor,
    lambda          = lambda_seq,
    kappa           = kappa_seq,
    ...
  )
  
  out <- list(
    nfolds    = nfolds,
    folds     = foldid,
    lambda    = lambda_seq,
    kappa     = kappa_seq,
    pairs     = pairs,
    fold_loss = fold_loss,
    cv_total  = cv_total,
    cv_mean   = cv_mean,
    idx_min   = idx_min,
    fit       = final_fit
  )
  class(out) <- "cv_glmtlp_ce"
  out
}

# S3 method so you can keep using coef(cv_obj)
coef.cv_glmtlp_ce <- function(object, which = object$idx_min, ...) {
  res <- tryCatch(
    predict(object$fit, type = "coefficients", which = which),
    error = function(e) {
      predict(object$fit, type = "coefficients", s = which)
    }
  )
  drop(res)
}

## ===============================================================
## 0 — working directory  (edit if needed)
## ===============================================================
#setwd("/Users/hongruzhao/myfile/Research simulation/Federated learning/data_final")
setwd("~/myfile/Research simulation/Federated learning/data_final")

## ===============================================================
## 1 — helper: load one hospital & append `site`
## ===============================================================
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

## ===============================================================
## 2 — load all four hospitals
## ===============================================================
sites <- c("cleveland", "hungarian", "swiss", "longbeach")
datasets <- lapply(setNames(sites, sites),
                   load_one_site,
                   all_sites = sites)

## quick check
sapply(datasets, function(d) dim(d$X))   # rows × 14 (13 + site)

## ===============================================================
## 3 — stack all X, y
## ===============================================================
X_all <- do.call(rbind, lapply(datasets, `[[`, "X"))
rownames(X_all) <- NULL
y_all <- unlist(lapply(datasets, `[[`, "y"), use.names = FALSE)

## ===============================================================
## 4 — split variable roles
## ===============================================================
global_vars             <- c("age","sex","chol","trestbps","thalach","oldpeak", "cp")
site_vars_lack_missing  <- c("fbs","restecg","exang")
site_vars_huge_missing  <- c("slope","ca","thal")
site_vars               <- c(site_vars_lack_missing, site_vars_huge_missing)

## keep global predictors ***plus site***
X_global <- X_all[ , c(global_vars, "site") ]

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

## ===============================================================
## 5 — separate Cleveland vs non-Cleveland
## ===============================================================
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

## ===============================================================
#################################################################
##  6 — build numeric design matrices with uniform column order
#################################################################

## --------------------------------------------------
## 6.1 GLOBAL block (same for every row)
## --------------------------------------------------
form_global <- ~   age + site + sex + chol + trestbps + thalach + oldpeak + cp
Z_global_all <- model.matrix(form_global, data = rbind(X_cle_noNA[,-13], X_noncle_noNA))
Z_global_all <- Z_global_all[,-1]  # drop intercept
colnames(Z_global_all)

## --------------------------------------------------
## 6.2 SITE-SPECIFIC block (site main-effects + interactions)
## --------------------------------------------------
form_site_cle <- ~ -1 + (fbs + restecg + exang + slope + ca + thal)
Z_site_cle <- model.matrix(form_site_cle, data = X_cle_noNA )

form_site_rest <- ~ -1 + site:(fbs + restecg + exang + slope + thal)
Z_site_rest <- model.matrix(form_site_rest, data = X_noncle_noNA )

## if a column name contains “cleveland”, remove it
clev_mask <- grepl("cleveland", colnames(Z_site_rest), fixed = TRUE)
if (any(clev_mask)) {
  Z_site_rest <- Z_site_rest[ , !clev_mask, drop = FALSE ]
}

## 7·1 split Z_global_all into the two row groups
n_cle   <- nrow(Z_site_cle)
n_rest  <- nrow(Z_site_rest)

Z_global_cle   <- Z_global_all[  1:n_cle,  , drop = FALSE ]
Z_global_rest  <- Z_global_all[(n_cle+1):nrow(Z_global_all), , drop = FALSE ]

## 7·2 placeholders so col-counts match
zero_cle_block  <- matrix(0, nrow = n_cle,  ncol = ncol(Z_site_rest),
                          dimnames = list(NULL, colnames(Z_site_rest)))
zero_rest_block <- matrix(0, nrow = n_rest, ncol = ncol(Z_site_cle),
                          dimnames = list(NULL, colnames(Z_site_cle)))

## 7·3 bind columns for each row group
Z_cle_final    <- cbind(Z_global_cle,  Z_site_cle,  zero_cle_block)
Z_rest_final   <- cbind(Z_global_rest, zero_rest_block, Z_site_rest)

## 7·4 stack into one grand design matrix
Z_all_final <- rbind(Z_cle_final, Z_rest_final)

y_final <- c(y_cle_noNA , y_noncle_noNA)

## ===============================================================
## Allowed tokens (for convenience in hetero_LRT)
## ===============================================================
ALLOWED_GLOBALS <- c("age","sex","chol","trestbps","thalach","oldpeak","cp","site")
ALLOWED_TOKENS  <- c(ALLOWED_GLOBALS,
                     # site-specific & interactions if ever needed in exploration:
                     "fbs","restecg","exang","slope","ca","thal",
                     "site:fbs","site:restecg","site:exang","site:slope","site:thal")

## ===============================================================
## hetero_LRT() — wrapper with clean output & print method
##   (keeps your current logic: CV only for H0; H1 = logistic GLM on union support)
## ===============================================================
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
  if (length(bad) > 0) {
    warning("Ignoring invalid globals: ", paste(bad, collapse = ", "),
            "\nAllowed: ", paste(allowed, collapse = ", "))
    globals <- setdiff(globals, bad)
    if (length(globals) == 0) stop("No valid globals remain after filtering.")
  }
  
  if (family != "binomial") stop("hetero_LRT currently supports family='binomial' only.")
  if (penalty != "l0") stop("hetero_LRT currently supports penalty='l0' only.")
  
  ## 2) Identify GLOBAL block columns to test (prefix match)
  pat <- paste0("^(", paste(globals, collapse = "|"), ")")
  B_index <- grep(pat, colnames(Z))
  if (length(B_index) == 0)
    stop("No matching columns in the design matrix for: ",
         paste(globals, collapse = ", "), call. = FALSE)
  
  ## 3) Build H0 / H1 design matrices
  Z_H1 <- Z
  Z_H0 <- Z[, -B_index, drop = FALSE]
  
  ## 4) Cross-validated fit (H0 ONLY) using cv_glmtlp_ce
  set.seed(seed)
  cv_H0 <- cv_glmtlp_ce(
    X = Z_H0, y = y,
    nfolds = nfolds, standardize = FALSE
  )
  idxH0 <- cv_H0$idx_min
  lamH0 <- cv_H0$pairs$lambda[idxH0]
  kapH0 <- cv_H0$pairs$kappa[idxH0]
  
  ## 5) Supports (from cv_H0) and union with B
  b0_cv <- coef(cv_H0)  # named vector: "(Intercept)", terms...
  nzcoef_H0 <- b0_cv[b0_cv != 0]
  support_H0 <- setdiff(names(nzcoef_H0), "(Intercept)")
  
  supp_H0_names <- support_H0
  supp_union_names <- union(supp_H0_names, colnames(Z_H1)[B_index])
  supp_union_idx   <- match(supp_union_names, colnames(Z_H1))
  supp_union_idx   <- sort(supp_union_idx[!is.na(supp_union_idx)])
  
  ## 6) Logistic GLM on union support (H1)
  X_H1_union <- if (length(supp_union_idx)) Z_H1[, supp_union_idx, drop = FALSE] else
    matrix(numeric(0), nrow = nrow(Z_H1), ncol = 0)
  
  fit_glm_H1 <- glm.fit(
    x       = cbind("(Intercept)" = 1, X_H1_union),
    y       = y,
    family  = binomial()
  )
  b1_glm <- fit_glm_H1$coefficients
  b1_glm[is.na(b1_glm)] <- 0
  nzcoef_H1 <- b1_glm[b1_glm != 0]
  support_H1 <- setdiff(names(nzcoef_H1), "(Intercept)")
  
  ## 7) Likelihood-ratio test using current logic
  log1pexp <- function(z) ifelse(z > 30, z, log1p(exp(z)))  # numerical guard
  
  eta0 <- as.vector(Z_H0 %*% b0_cv[-1] + b0_cv[1])          # H0: cv-ℓ0 coefficients
  eta1 <- as.vector(cbind(1, X_H1_union) %*% b1_glm)        # H1: logistic GLM on union support
  
  ll0 <- sum(y * eta0 - log1pexp(eta0))
  ll1 <- sum(y * eta1 - log1pexp(eta1))
  
  LR_stat <- -2 * (ll0 - ll1)
  df      <- length(B_index)
  p_val   <- pchisq(LR_stat, df = df, lower.tail = FALSE)
  
  ## 8) Return structured result (with names & supports)
  out <- list(
    tested_globals     = globals,
    df                 = df,
    LR                 = as.numeric(LR_stat),
    p_value            = as.numeric(p_val),
    B_index            = B_index,
    B_columns          = colnames(Z)[B_index],
    # CV-H0 tuning snapshot
    selected_idx_H0    = idxH0,
    selected_lambda_H0 = lamH0,
    selected_kappa_H0  = kapH0,
    # objects & coefficients
    cv_H0              = cv_H0,
    coef_H0            = b0_cv,
    coef_H1            = b1_glm,
    # supports & nonzero coefs
    nzcoef_H0          = nzcoef_H0,
    nzcoef_H1          = nzcoef_H1,
    support_H0         = support_H0,
    support_H1         = support_H1,
    support_union      = supp_union_names
  )
  class(out) <- c("hetero_LRT_result", class(out))
  out
}

## Pretty print
print.hetero_LRT_result <- function(x, max_terms = 30, ...) {
  cat("Heterogeneous LR test (binomial; H0 via CV-L0, H1 via GLM on union support)\n")
  cat("  Tested globals :", paste(x$tested_globals, collapse = ", "), "\n")
  cat("  df             :", x$df, "\n")
  cat("  LR statistic   :", signif(x$LR, 5), "\n")
  cat("  p-value        :", signif(x$p_value, 5), "\n")
  if (!is.null(x$B_columns)) {
    cat("  Columns tested :", paste(x$B_columns, collapse = ", "), "\n")
  }
  cat("  H0 CV choice   : idx =", x$selected_idx_H0,
      ", lambda =", signif(x$selected_lambda_H0, 4),
      ", kappa =", signif(x$selected_kappa_H0, 4), "\n")
  
  prettify <- function(v, k = max_terms) {
    if (length(v) == 0) return("<none>")
    if (length(v) <= k)  return(paste(v, collapse = ", "))
    paste(paste(v[seq_len(k)], collapse = ", "),
          sprintf("... (+%d more)", length(v) - k))
  }
  
  cat("\n  Support(H0) size:", length(x$support_H0), "\n")
  cat("  Support(H0)     :", prettify(x$support_H0), "\n")
  cat("  Support(H1) size:", length(x$support_H1)+1, "\n")
  cat("  Support(H1)     :", prettify(x$support_H1), "\n")
  cat("  Union support   :", prettify(x$support_union), "\n")
  
  cat("\n  Nonzero coefficients under H0 (incl. intercept):\n")
  print(x$nzcoef_H0)
  cat("\n  Nonzero coefficients under H1 (incl. intercept):\n")
  print(x$nzcoef_H1)
  
  invisible(x)
}

## ===============================================================
## 9 — Minimal driver (examples) — uncomment to run
## ===============================================================
# Example single tests
res_age  <- hetero_LRT("age")
res_sex  <- hetero_LRT("sex")
res_chol <- hetero_LRT("chol")
res_cp   <- hetero_LRT("cp")
print(res_age); 
print(res_sex); 
print(res_chol); 
print(res_cp)

# Example: multiple globals jointly
res_two <- hetero_LRT(c("age","sex"))
print(res_two)

# Batch over allowed globals (collect a small summary table)
results_all <- lapply(ALLOWED_GLOBALS, function(v) {
   r <- hetero_LRT(v)
   data.frame(var = v, df = r$df, LR = r$LR, p = r$p_value)
 })
do.call(rbind, results_all)

