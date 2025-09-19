## ──────────────────────────────────────────────────────────────────────────────
##  0. Libraries & seed
## ──────────────────────────────────────────────────────────────────────────────
library(MASS)      # mvrnorm
library(glmtlp)    # L0-penalised GLM
set.seed(2026)

## ──────────────────────────────────────────────────────────────────────────────
##  1. Helper functions
## ──────────────────────────────────────────────────────────────────────────────

## 1.1  Design generator (equal-correlation, separate X and W)
generate_design <- function(K, n, p0, pj, rho) {
  max_p <- max(p0, pj)
  Sigma <- matrix(rho, max_p, max_p); diag(Sigma) <- 1
  
  lapply(seq_len(K), function(j) {
    ## global block X: intercept + MVN(p0−1)
    X <- if (p0 > 1) {
      X_rand <- MASS::mvrnorm(n, rep(0, p0 - 1),
                              Sigma[1:(p0 - 1), 1:(p0 - 1)])
      cbind(Intercept = 1, X_rand)
    } else {
      matrix(1, n, 1, dimnames = list(NULL, "Intercept"))
    }
    ## site-specific block W
    W <- MASS::mvrnorm(n, rep(0, pj), Sigma[1:pj, 1:pj])
    list(X = X, W = W)
  })
}

## 1.2  β generators
generate_beta0 <- function(p0, B_idx, delta) {
  beta0 <- numeric(p0); beta0[1:3] <- c(1, 2, 3)
  beta0[B_idx[1]] <- delta
  beta0
}
generate_beta_list <- function(p_vec) {
  lapply(seq_along(p_vec), function(j) {
    b <- numeric(p_vec[j]); b[1] <- j; b
  })
}

## 1.3  Data generator
generate_data <- function(design, beta0, beta_list, sigma_error = 1) {
  mapply(function(site, b_j) {
    eps <- rnorm(nrow(site$X), sd = sigma_error)
    as.vector(site$X %*% beta0 + site$W %*% b_j + eps)
  }, design, beta_list, SIMPLIFY = FALSE)
}

## 1.4  Assemble full block-diagonal design
assemble_full_model <- function(design, beta0, beta_list, Y_list) {
  W_list <- lapply(design, `[[`, "W");   p_vec <- sapply(W_list, ncol)
  row_blocks <- lapply(seq_along(design), function(j) {
    Xj <- design[[j]]$X;  Wj <- W_list[[j]];  nj <- nrow(Xj)
    left  <- if (j > 1) matrix(0, nj, sum(p_vec[1:(j - 1)])) else NULL
    right <- if (j < length(design)) matrix(0, nj, sum(p_vec[(j + 1):length(design)])) else NULL
    cbind(Xj, left, Wj, right)
  })
  list(
    design   = do.call(rbind, row_blocks),
    beta     = c(beta0, unlist(beta_list)),
    response = unlist(Y_list)
  )
}

## ──────────────────────────────────────────────────────────────────────────────
##  2. Simulation settings
## ──────────────────────────────────────────────────────────────────────────────
p_grid    <- 100                  # total p
B_sizes   <- 1
delta_vec <- c(0, 0.1, 0.2)       # signal sizes
n_vec     <- rep(200, 3)          # per-site training size
p_vec     <- c(10, 10, 10)        # per-site W dims
rho_eq    <- 0.05                 # equal correlation
n_rep     <- 400                  # Monte-Carlo replications
K         <- length(n_vec)

## validation pool (independent)
n_val_total <- 3000
n_val_site  <- ceiling(n_val_total / K)

results <- list()

## ──────────────────────────────────────────────────────────────────────────────
##  3. Main simulation loop
## ──────────────────────────────────────────────────────────────────────────────
p0 <- p_grid - sum(p_vec)

for (B_size in B_sizes) {
  B_idx <- 4:(3 + B_size)                         # indices before dropping intercept
  
  for (delta in delta_vec) {
    rejections <- logical(n_rep)
    
    for (rep in seq_len(n_rep)) {
      ## 3.1  Generate *training* data
      tr_design <- lapply(seq_len(K), function(j)
        generate_design(1, n_vec[j], p0, p_vec[j], rho_eq)[[1]])
      beta0_true <- generate_beta0(p0, B_idx, delta)
      beta_list  <- generate_beta_list(p_vec)
      Y_tr_list  <- generate_data(tr_design, beta0_true, beta_list, sigma_error = 0.1)
      tr_full    <- assemble_full_model(tr_design, beta0_true, beta_list, Y_tr_list)
      X_tr       <- tr_full$design
      Y_tr       <- tr_full$response
      n_tr       <- length(Y_tr)
      
      ## 3.2  Generate *validation* data
      val_design <- lapply(seq_len(K), function(j)
        generate_design(1, n_val_site, p0, p_vec[j], rho_eq)[[1]])
      Y_val_list <- generate_data(val_design, beta0_true, beta_list, sigma_error = 0.1)
      val_full   <- assemble_full_model(val_design, beta0_true, beta_list, Y_val_list)
      X_val      <- val_full$design
      Y_val      <- val_full$response
      
      ## 3.3  Drop intercept column (glmtlp will estimate its own)
      X_tr2  <- X_tr [, -1, drop = FALSE]
      X_val2 <- X_val[, -1, drop = FALSE]
      B_shift <- B_idx - 1    # adjust indices after intercept removal
      
      ## ─── H0  : exclude B columns, fit λ-path, tune on validation ───────────
      Xtr_H0  <- X_tr2 [, -B_shift, drop = FALSE]
      Xval_H0 <- X_val2[, -B_shift, drop = FALSE]
      
      fit_H0 <- glmtlp(
        Xtr_H0, Y_tr,
        family  = "gaussian",
        penalty = "l0",
        nlambda = 100,
        lambda.min.ratio = 1e-3
      )
      preds_H0 <- predict(fit_H0, Xval_H0)           # n_val × nlambda matrix
      mse_H0   <- colMeans((Y_val - preds_H0)^2)
      best_H0  <- which.min(mse_H0)
      yhat_H0  <- as.vector(predict(fit_H0, Xtr_H0, s = fit_H0$lambda[best_H0]))
      
      ## ─── H1  : full design, un-penalise B columns, tune on validation ──────
      pfac          <- rep(1, ncol(X_tr2))
      pfac[B_shift] <- 0
      fit_H1 <- glmtlp(
        X_tr2, Y_tr,
        family  = "gaussian",
        penalty = "l0",
        pfactor = pfac,
        nlambda = 100,
        lambda.min.ratio = 1e-3
      )
      preds_H1 <- predict(fit_H1, X_val2)
      mse_H1   <- colMeans((Y_val - preds_H1)^2)
      best_H1  <- which.min(mse_H1)
      yhat_H1  <- as.vector(predict(fit_H1, X_tr2, s = fit_H1$lambda[best_H1]))
      
      ## 3.4  Compute LRT
      rss0      <- sum((Y_tr - yhat_H0)^2)
      rss1      <- sum((Y_tr - yhat_H1)^2)
      sigma0_sq <- rss0 / n_tr
      sigma1_sq <- rss1 / n_tr
      LRT       <- n_tr * log(sigma0_sq / sigma1_sq)
      
      rejections[rep] <- (LRT > qchisq(0.95, df = length(B_idx)))
    } # rep
    
    key <- paste(p_grid, B_size, delta, sep = "_")
    results[[key]] <- mean(rejections)
  } # delta
}   # B_size

## ──────────────────────────────────────────────────────────────────────────────
##  4. Post-processing
## ──────────────────────────────────────────────────────────────────────────────
results_df <- do.call(rbind, lapply(names(results), function(k) {
  parts <- as.numeric(strsplit(k, "_")[[1]])
  data.frame(
    p              = parts[1],
    B_size         = parts[2],
    delta          = parts[3],
    rejection_rate = results[[k]]
  )
}))
rownames(results_df) <- NULL
print(results_df)
