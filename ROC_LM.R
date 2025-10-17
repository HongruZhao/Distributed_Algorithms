## --- 0.  Set-up -------------------------------------------------------

library(MASS)
library(glmtlp)
set.seed(2025)

## --- 1.  Data-generation helpers -------------------------------------

### 1.1  Equal-correlation design (col-1 = 1)
generate_design <- function(K, n, p0, pj, rho) {
  max_p <- max(p0, pj)
  Σ <- matrix(rho, max_p, max_p); diag(Σ) <- 1
  
  lapply(seq_len(K), function(j) {
    X <- if (p0 > 1) cbind(1, mvrnorm(n, rep(0, p0 - 1), Σ[1:(p0 - 1), 1:(p0 - 1)]))
    else       matrix(1, n, 1)
    W <- mvrnorm(n, rep(0, pj), Σ[1:pj, 1:pj])
    list(X = X, W = W)
  })
}

### 1.2  Sparse β
generate_betas <- function(K, p0, pj, s0, sj, βmin, βmax = 2 * βmin) {
  β0 <- numeric(p0); β0[sample(p0, s0)] <- runif(s0, βmin, βmax) * sample(c(-1, 1), s0, TRUE)
  βlist <- lapply(seq_len(K), function(j) {
    b <- numeric(pj); b[sample(pj, sj)] <- runif(sj, βmin, βmax) * sample(c(-1, 1), sj, TRUE); b
  })
  list(beta0 = β0, beta_list = βlist)
}

### 1.3  Gaussian responses
generate_data <- function(designs, β0, βlist, σ2) {
  K <- length(designs); σ2 <- if (length(σ2) == 1) rep(σ2, K) else σ2
  lapply(seq_len(K), function(j) {
    X <- designs[[j]]$X; W <- designs[[j]]$W
    y <- X %*% β0 + W %*% βlist[[j]] + rnorm(nrow(X), 0, sqrt(σ2[j]))
    list(X = X, W = W, Y = as.numeric(y))
  })
}

### 1.4  Assemble block-matrix   (Intercept col kept, later dropped)
assemble_block_data <- function(sim_data, betas) {
  K   <- length(sim_data)
  p0  <- length(betas$beta0); pj <- length(betas$beta_list[[1]])
  n_j <- sapply(sim_data, \(d) nrow(d$X))
  
  N <- sum(n_j); P <- p0 + K * pj
  Xall <- matrix(0, N, P); Yall <- numeric(N)
  
  start <- 1
  for (j in seq_len(K)) {
    rows <- start:(start + n_j[j] - 1)
    Xall[rows, 1:p0] <- sim_data[[j]]$X
    cs <- p0 + (j - 1) * pj + 1; ce <- p0 + j * pj
    Xall[rows, cs:ce] <- sim_data[[j]]$W
    Yall[rows] <- sim_data[[j]]$Y
    start <- start + n_j[j]
  }
  list(X = Xall, Y = Yall, beta = c(betas$beta0, unlist(betas$beta_list)))
}

## --- 2.  Main simulation+estimation function --------------------------

run_simulation_cv <- function(params) {
  with(params, {
    roc_list <- vector("list", n_reps)
    
    for (r in seq_len(n_reps)) {
      ## 2.1  Generate one replication
      des  <- generate_design(K, n_val, p0_val, pj_val, rho_val)
      bts  <- generate_betas(K, p0_val, pj_val, s0_val, sj_val, beta_min_val)
      simd <- generate_data(des, bts$beta0, bts$beta_list, sigma2_val)
      bd   <- assemble_block_data(simd, bts)
      
      ## 2.2  Drop the hand-made intercept column
      X <- bd$X[, -1, drop = FALSE]; Y <- bd$Y
      N <- nrow(X)
      
      ## 2.3  50/50 split
      tr_idx <- sample(N, N/2); va_idx <- setdiff(seq_len(N), tr_idx)
      X_tr <- X[tr_idx, , drop = FALSE]; Y_tr <- Y[tr_idx]
      X_va <- X[va_idx, , drop = FALSE]; Y_va <- Y[va_idx]
      
      ## 2.4  Fit L0 path on train
      fit <- glmtlp(
        X_tr, Y_tr,
        family = "gaussian", penalty = "l0",
        nlambda = nlambda, lambda.min.ratio = 1e-3,
        intercept = TRUE
      )
      Kgrid <- length(fit$kappa)
      
      ## 2.5  Choose λ with min validation-MSE
      Yhat_va <- predict(fit, X = X_va, type = "response")
      best_k  <- which.min(colMeans((Y_va - Yhat_va)^2))
      
      ## 2.6  Coefficient path
      coef_mat <- coef(fit, which = seq_len(Kgrid), drop = FALSE)    # (1+P)×Kgrid
      coef_pred <- coef_mat[-1, , drop = FALSE]                      # drop intercept
      beta_star <- coef_pred[, best_k]
      
      ## 2.7  True support (predictors only)
      true_beta <- bd$beta[-1]
      tsupp     <- which(true_beta != 0)
      nulls     <- setdiff(seq_len(nrow(coef_pred)), tsupp)
      
      ## 2.8  Full ROC path
      TPR <- FPR <- numeric(Kgrid)
      for (k in seq_len(Kgrid)) {
        sel <- which(abs(coef_pred[, k]) > 1e-8)
        TPR[k] <- length(intersect(sel, tsupp)) / length(tsupp)
        FPR[k] <- length(setdiff(sel, tsupp))   / length(nulls)
      }
      
      ## 2.9  Tuned point
      sel_star <- which(abs(beta_star) > 1e-8)
      TPR_star <- length(intersect(sel_star, tsupp)) / length(tsupp)
      FPR_star <- length(setdiff(sel_star, tsupp))   / length(nulls)
      
      roc_list[[r]] <- list(FPR = FPR, TPR = TPR,
                            FPR_star = FPR_star, TPR_star = TPR_star)
    }
    roc_list
  })
}

## --- 3.  Parameter blocks --------------------------------------------

params1 <- list(
  K = 3, n_val = 200, p0_val = 50, pj_val = 50,
  s0_val = 10, sj_val = 10, rho_val = 0.6,
  beta_min_val = 0.4, sigma2_val = rep(2, 3),
  n_reps = 10, nlambda = 100
)

params2 <- list(
  K = 3, n_val = 200, p0_val = 50, pj_val = 50,
  s0_val = 10, sj_val = 10, rho_val = 0.3,
  beta_min_val = 0.6, sigma2_val = rep(1, 3),
  n_reps = 10, nlambda = 100
)

## --- 4.  Run simulations ---------------------------------------------

roc1 <- run_simulation_cv(params1)
roc2 <- run_simulation_cv(params2)

## --- 5.  Plotting -----------------------------------------------------

pdf("roc_plots_LM.pdf", 10, 5)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

## Panel 1
plot(roc1[[1]]$FPR, roc1[[1]]$TPR, type = "l", col = "grey",
     xlab = "False Positive Rate", ylab = "True Positive Rate")
for (r in 2:length(roc1)) lines(roc1[[r]]$FPR, roc1[[r]]$TPR, col = "grey")
for (r in seq_along(roc1))
  points(roc1[[r]]$FPR_star, roc1[[r]]$TPR_star, pch = 4, col = "blue")
legend("bottomright", legend = "Validation-chosen", pch = 4,
       col = "blue", bty = "n")

## Panel 2
plot(roc2[[1]]$FPR, roc2[[1]]$TPR, type = "l", col = "grey",
     xlab = "False Positive Rate", ylab = "True Positive Rate")
for (r in 2:length(roc2)) lines(roc2[[r]]$FPR, roc2[[r]]$TPR, col = "grey")
for (r in seq_along(roc2))
  points(roc2[[r]]$FPR_star, roc2[[r]]$TPR_star, pch = 4, col = "blue")
legend("bottomright", legend = "Validation-chosen", pch = 4,
       col = "blue", bty = "n")

dev.off()

