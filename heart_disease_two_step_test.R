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

  # glmtlp objects (especially for L0) report a 2-D grid; the number of
  # solutions we can actually evaluate is the length of the intercept path.
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
        # fallback to common 'newx'/'s' arguments if needed
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
      # ensure columns correspond to indices 1..L
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

## grab individual objects if you want them
X_cleveland <- datasets$cleveland$X; y_cleveland <- datasets$cleveland$y
X_hungarian <- datasets$hungarian$X; y_hungarian <- datasets$hungarian$y
X_swiss     <- datasets$swiss$X;     y_swiss     <- datasets$swiss$y
X_longbeach <- datasets$longbeach$X; y_longbeach <- datasets$longbeach$y

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
X_site_huge_missing_non_cleveland[] <-
  lapply(X_site_huge_missing_non_cleveland, factor)

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
form_global <- ~   age +site+ sex + chol + trestbps + thalach + oldpeak+cp
Z_global_all <- model.matrix(form_global, data = rbind(X_cle_noNA[,-13], X_noncle_noNA))
Z_global_all <- Z_global_all[,-1]
colnames(Z_global_all)

## --------------------------------------------------
## 6.2 SITE-SPECIFIC block (site main-effects + interactions)
## --------------------------------------------------
form_site_cle <- ~-1 + (  fbs + restecg + exang + slope + ca + thal)
Z_site_cle <- model.matrix(form_site_cle, data = X_cle_noNA )
colnames(Z_site_cle)

form_site_rest <- ~-1 + site:(  fbs + restecg + exang + slope + thal)
Z_site_rest <- model.matrix(form_site_rest, data = X_noncle_noNA )
colnames(Z_site_rest)

## if a column name contains “cleveland”, remove it
clev_mask <- grepl("cleveland", colnames(Z_site_rest), fixed = TRUE)
if (any(clev_mask)) {
  Z_site_rest <- Z_site_rest[ , !clev_mask, drop = FALSE ]
}
colnames(Z_site_rest)         # should show no “cleveland” columns

## 7·1  split Z_global_all into the two row groups
n_cle   <- nrow(Z_site_cle)
n_rest  <- nrow(Z_site_rest)

Z_global_cle   <- Z_global_all[  1:n_cle,  , drop = FALSE ]
Z_global_rest  <- Z_global_all[(n_cle+1):nrow(Z_global_all), , drop = FALSE ]

## 7·2  build zero-filled placeholders so col-counts match
zero_cle_block  <- matrix(0, nrow = n_cle,  ncol = ncol(Z_site_rest),
                          dimnames = list(NULL, colnames(Z_site_rest)))
zero_rest_block <- matrix(0, nrow = n_rest, ncol = ncol(Z_site_cle),
                          dimnames = list(NULL, colnames(Z_site_cle)))

## 7·3  bind columns for each row group
Z_cle_final    <- cbind(Z_global_cle,  Z_site_cle,  zero_cle_block)
Z_rest_final   <- cbind(Z_global_rest, zero_rest_block, Z_site_rest)

## 7·4  stack into one grand design matrix
Z_all_final <- rbind(Z_cle_final, Z_rest_final)

y_final <- c(y_cle_noNA , y_noncle_noNA)

## ===============================================================
## 8 — Specify B-set and H0/H1 design matrices
## ===============================================================
library(glmtlp)

# Example: test only the GLOBAL variable 'age' (edit if needed)
# global_vars <- c("age","sex","chol","trestbps","thalach","oldpeak")
global_vars <- c("age")

B_index <- grep(paste0("^(", paste(global_vars, collapse="|"), ")"),
                colnames(Z_all_final))
cat("Testing indices B_index =", B_index, "\n",
    "columns:\n", colnames(Z_all_final)[B_index], "\n\n")

# H1: full heterogeneous model
Z_H1 <- Z_all_final

# H0: GLOBAL block removed → β_B = 0
Z_H0 <- Z_all_final[, -B_index, drop = FALSE]

## ===============================================================
## 9 — 5-fold CV fit under H0 and H1 (REPLACING cv.glmtlp)
## ===============================================================
set.seed(2025)

# (Optional) A single overall CV just to inspect the path on the full design
cv_model <- cv_glmtlp_ce(
  X = Z_all_final, y = y_final,
  nfolds = 5, standardize = FALSE
)
cat("Selected idx (all):", cv_model$idx_min,
    " lambda:", cv_model$pairs$lambda[cv_model$idx_min],
    " kappa:",  cv_model$pairs$kappa[cv_model$idx_min], "\n")
# Nonzero coefficients at the selected index:
coef(cv_model)[coef(cv_model) != 0]

# H0 fit: all columns penalized (default)
cv_H0 <- cv_glmtlp_ce(
  X = Z_H0, y = y_final,
  nfolds = 5, standardize = FALSE
)
cat("Selected idx (H0):", cv_H0$idx_min,
    " lambda:", cv_H0$pairs$lambda[cv_H0$idx_min],
    " kappa:",  cv_H0$pairs$kappa[cv_H0$idx_min], "\n")
coef(cv_H0)[coef(cv_H0) != 0]

# H1 fit: unpenalize the GLOBAL block only here
p   <- ncol(Z_H1)
pf  <- rep(1, p)
pf[B_index] <- 0

cv_H1 <- cv_glmtlp_ce(
  X = Z_H1, y = y_final,
  nfolds = 5, standardize = FALSE,
  penalty.factor = pf
)
cat("Selected idx (H1):", cv_H1$idx_min,
    " lambda:", cv_H1$pairs$lambda[cv_H1$idx_min],
    " kappa:",  cv_H1$pairs$kappa[cv_H1$idx_min], "\n")
coef(cv_H1)[coef(cv_H1) != 0]

## ===============================================================
## 10 — Likelihood-ratio test with selected fits
## ===============================================================
# Coef vectors (Intercept first, then columns of the respective X)
b0 <- coef(cv_H0)
b1 <- coef(cv_H1)

# Linear predictors
eta0 <- as.vector(Z_H0 %*% b0[-1] + b0[1])
eta1 <- as.vector(Z_H1 %*% b1[-1] + b1[1])

# Log-likelihoods ℓ(β) = ∑[y η − log(1+e^η)]
ll0 <- sum(y_final * eta0 - log1p(exp(eta0)))
ll1 <- sum(y_final * eta1 - log1p(exp(eta1)))

# LR statistic and chi-square p-value (df = |B_index|)
LR_stat <- -2 * (ll0 - ll1)
df      <- length(B_index)
p_val   <- pchisq(LR_stat, df = df, lower.tail = FALSE)
cat("LR =", round(LR_stat, 3), " df =", df, " p =", signif(p_val, 3), "\n")
