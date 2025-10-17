## --- 0.  Working directory & packages --------------------------------

# ----------------------------------------------------------
# Data generation • 50/50 split • L0-logistic • CE tuning • ROC plotting
# ----------------------------------------------------------



library(MASS)
library(glmtlp)
set.seed(2025)                             # reproducibility

## --- 1.  Helper functions ---------------------------------------------

### 1.1  Equal-correlation design (first X column = 1)
generate_design <- function(K, n, p0, pj, rho) {
  max_p <- max(p0, pj)
  Sigma <- matrix(rho, max_p, max_p); diag(Sigma) <- 1
  
  lapply(seq_len(K), function(j) {
    X <- if (p0 > 1) {
      cbind(1, mvrnorm(n, rep(0, p0 - 1), Sigma[1:(p0 - 1), 1:(p0 - 1)]))
    } else {
      matrix(1, n, 1)
    }
    W <- mvrnorm(n, rep(0, pj), Sigma[1:pj, 1:pj])
    list(X = X, W = W)
  })
}

### 1.2  Random β with sparsity
generate_betas <- function(K, p0, pj, s0, sj, beta_min, beta_max = 2 * beta_min) {
  beta0 <- numeric(p0)
  beta0[sample(p0, s0)] <- runif(s0, beta_min, beta_max) * sample(c(-1, 1), s0, TRUE)
  beta_list <- lapply(seq_len(K), function(j) {
    bj <- numeric(pj)
    bj[sample(pj, sj)] <- runif(sj, beta_min, beta_max) * sample(c(-1, 1), sj, TRUE)
    bj
  })
  list(beta0 = beta0, beta_list = beta_list)
}

### 1.3  Simulate Y ~ Bernoulli(logit⁻¹(Xβ + Wβ_j))
generate_data <- function(design_list, beta0, beta_list) {
  lapply(seq_along(design_list), function(j) {
    X <- design_list[[j]]$X
    W <- design_list[[j]]$W
    p <- plogis(X %*% beta0 + W %*% beta_list[[j]])
    list(X = X, W = W, Y = rbinom(nrow(X), 1, p))
  })
}

### 1.4  Assemble block matrix
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

## --- 2.  50/50-split, L0-logistic, CE-tuned ROC -----------------------

run_simulation_cv_logistic <- function(params) {
  with(params, {
    roc_list <- vector("list", n_reps)
    
    for (r in seq_len(n_reps)) {
      ## 2.1  Simulate one replication
      des  <- generate_design(K, n_val, p0_val, pj_val, rho_val)
      bts  <- generate_betas(K, p0_val, pj_val, s0_val, sj_val, beta_min_val)
      simd <- generate_data(des, bts$beta0, bts$beta_list)
      bd   <- assemble_block_data(simd, bts)
      
      ## 2.2  Drop the manual intercept column
      X <- bd$X[, -1, drop = FALSE]; Y <- bd$Y
      N <- nrow(X)
      
      ## 2.3  60/40 split
      tr_idx <- sample(N, size = floor(0.5 * N))
      va_idx <- setdiff(seq_len(N), tr_idx)
      X_tr <- X[tr_idx, , drop = FALSE]; Y_tr <- Y[tr_idx]
      X_va <- X[va_idx, , drop = FALSE]; Y_va <- Y[va_idx]
      
      ## 2.4  Fit λ-path on training
      fit <- glmtlp(
        X_tr, Y_tr,
        family           = "binomial",
        penalty          = "l0",
        nlambda          = nlambda,
        lambda.min.ratio = 1e-3,
        intercept        = TRUE
      )
      Kgrid <- length(fit$kappa)
      
      ## 2.5  Pick λ by minimum CE on validation
      P_hat <- predict(fit, X = X_va, type = "response")
      P_hat <- pmin(pmax(P_hat, 1e-8), 1 - 1e-8)     # clip
      CE    <- colMeans(-(Y_va * log(P_hat) + (1 - Y_va) * log(1 - P_hat)))
      best_k <- which.min(CE)
      
      ## 2.6  Coefficients & true support
      coefs_full <- coef(fit, which = seq_len(Kgrid), drop = FALSE)
      coefs_n    <- coefs_full[-1, , drop = FALSE]       # drop intercept
      best_beta  <- coefs_n[, best_k]
      
      true_beta  <- bd$beta[-1]
      tsupp      <- which(true_beta != 0)
      nulls      <- setdiff(seq_len(nrow(coefs_n)), tsupp)
      
      ## 2.7  Full ROC path
      TPR <- FPR <- numeric(Kgrid)
      for (k in seq_len(Kgrid)) {
        sel <- which(abs(coefs_n[, k]) > 1e-8)
        TPR[k] <- length(intersect(sel, tsupp)) / length(tsupp)
        FPR[k] <- length(setdiff(sel, tsupp))   / length(nulls)
      }
      
      ## 2.8  Tuned (FPR*, TPR*)
      sel_star <- which(abs(best_beta) > 1e-8)
      TPR_star <- length(intersect(sel_star, tsupp)) / length(tsupp)
      FPR_star <- length(setdiff(sel_star, tsupp))   / length(nulls)
      
      ## 2.9  Store
      roc_list[[r]] <- list(
        FPR      = FPR,
        TPR      = TPR,
        FPR_star = FPR_star,
        TPR_star = TPR_star
      )
    }
    roc_list
  })
}

## --- 3.  Scenario parameters -----------------------------------------

params1 <- list(
  K = 3, n_val = 500, p0_val = 30, pj_val = 30,
  s0_val = 10, sj_val = 5, rho_val = 0.6,
  beta_min_val = 0.4, n_reps = 10, nlambda = 60
)

params2 <- list(
  K = 3, n_val = 500, p0_val = 30, pj_val = 30,
  s0_val = 10, sj_val = 5, rho_val = 0.2,
  beta_min_val = 1, n_reps = 10, nlambda = 60
)

## --- 4.  Run simulations ----------------------------------------------

roc1 <- run_simulation_cv_logistic(params1)
roc2 <- run_simulation_cv_logistic(params2)


## --- 5.  Plotting (skip any (0,0) points) -----------------------------
pdf("roc_logistic_ce.pdf", width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

draw_panel <- function(roc_list) {
  ## empty plot with fixed axes
  plot(NA, xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate", ylab = "True Positive Rate")
  
  ## grey ROC curves
  for (r in seq_along(roc_list)) {
    keep <- !(roc_list[[r]]$FPR == 0 & roc_list[[r]]$TPR == 0)
    if (any(keep))
      lines(roc_list[[r]]$FPR[keep], roc_list[[r]]$TPR[keep], col = "grey")
  }
  
  ## blue “×” for the tuned point (skip if it is (0,0))
  for (r in seq_along(roc_list)) {
    if (!(roc_list[[r]]$FPR_star == 0 && roc_list[[r]]$TPR_star == 0)) {
      points(roc_list[[r]]$FPR_star, roc_list[[r]]$TPR_star,
             pch = 4, col = "blue")
    }
  }
  
  legend("bottomright", legend = "Validation-chosen",
         pch = 4, col = "blue", bty = "n")
}

## Panel 1  (params 1)
draw_panel(roc1)

## Panel 2  (params 2)
draw_panel(roc2)

dev.off()
