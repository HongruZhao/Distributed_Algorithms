library(MASS)
library(glmtlp)
library(glmnet)

setwd("C:/HZ simulation")

set.seed(2025)

## ── 1. Helper functions ─────────────────────────────────────────────────────
generate_design <- function(K, n, p0, pj, rho) {
  max_p <- max(p0, pj)
  Sigma  <- matrix(rho, max_p, max_p); diag(Sigma) <- 1
  lapply(seq_len(K), function(j) {
    X <- if (p0 > 1) {
      cbind(Intercept = 1,
            MASS::mvrnorm(n, rep(0, p0 - 1),
                          Sigma[1:(p0 - 1), 1:(p0 - 1)]))
    } else matrix(1, n, 1, dimnames = list(NULL, "Intercept"))
    W <- MASS::mvrnorm(n, rep(0, pj), Sigma[1:pj, 1:pj])
    list(X = X, W = W)
  })
}

generate_beta0 <- function(p0, B_idx, delta, length_b0 = 3) {
  b0 <- numeric(p0)
  b0[1:length_b0] <- seq_len(length_b0)          # (1,2,3)
  b0[B_idx[1]]   <- delta                        # signal in *first* B entry
  b0
}

generate_beta_list  <- function(p_vec) lapply(seq_along(p_vec), \(j) {
  b <- numeric(p_vec[j]); b[1] <- j; b
})

generate_data <- function(design, beta0, beta_list, sigma = 1) {
  mapply(function(site, bj) {
    eps <- rnorm(nrow(site$X), sd = sigma)
    drop(site$X %*% beta0 + site$W %*% bj + eps)
  }, design, beta_list, SIMPLIFY = FALSE)
}

assemble_full_model <- function(design, beta0, beta_list, Y_list) {
  W_list <- lapply(design, `[[`, "W")
  p_vec  <- sapply(W_list, ncol)
  rows   <- lapply(seq_along(design), function(j) {
    Xj <- design[[j]]$X; Wj <- W_list[[j]]; nj <- nrow(Xj)
    left  <- if (j > 1) matrix(0, nj, sum(p_vec[1:(j - 1)])) else NULL
    right <- if (j < length(design)) matrix(0, nj, sum(p_vec[(j + 1):length(design)])) else NULL
    cbind(Xj, left, Wj, right)
  })
  list(design   = do.call(rbind, rows),
       beta     = c(beta0, unlist(beta_list)),
       response = unlist(Y_list))
}

## ── 2. Global settings ──────────────────────────────────────────────────────
p_grid       <- c(50, 100, 500, 1000)
B_sizes      <- c(1, 5, 10)
delta_vec    <- c(0.0, 0.1, 0.2, 0.3)
n_vec        <- rep(100, 3)
p_vec        <- c(10, 10, 10)
rho_eq       <- 0.5
n_rep        <- 1000
sigma_error  <- 1
length_b0    <- 3
true_sparse  <- length_b0 + 3 - 1
#kappa_0 <- as.integer(rep(true_sparse, 30))
kappa_0      <- rep(1:10, each = 3)   # (your updated grid)
K            <- length(n_vec)
alpha        <- 0.05                  # <<< CHANGED: set test level once

## ── 3. Containers for results ───────────────────────────────────────────────
results         <- list()
results_normal  <- list()
results_df_all  <- list()
row_id          <- 1L
results_banchmark <- list()  # will hold mean(rejections_banchmark)

## ── 4. Main simulation loops ────────────────────────────────────────────────
for (p in p_grid) {
  p0 <- p - sum(p_vec)
  for (B_size in B_sizes) {
    B_idx   <- (length_b0 + 1):(length_b0 + B_size)
    B_shift <- B_idx - 1
    
    for (delta in delta_vec) {
      
      rejections <- logical(n_rep)
      rejections_normal <- logical(n_rep)
      ## NEW: per-rep storage for the lasso-split LR benchmark (0/1)
      rejections_banchmark <- logical(n_rep)  # <<< CHANGED: keep only indicator
      
      for (rep in seq_len(n_rep)) {
        
        ## 4·1  Data generation ------------------------------------------------
        design_list <- lapply(seq_len(K), \(j)
                              generate_design(1, n_vec[j], p0, p_vec[j], rho_eq)[[1]])
        beta0_true <- generate_beta0(p0, B_idx, delta, length_b0)
        beta_list  <- generate_beta_list(p_vec)
        Y_list     <- generate_data(design_list, beta0_true, beta_list, sigma_error)
        
        full      <- assemble_full_model(design_list, beta0_true, beta_list, Y_list)
        X_full    <- full$design
        Y_full    <- full$response
        n_total   <- length(Y_full)
        X_noInt   <- X_full[, -1, drop = FALSE]
        
        ## 4·2·A  Sample-split LR benchmark (lasso on I1 excluding B; LR on I2)
        site_offsets <- c(0L, cumsum(n_vec))
        I1 <- unlist(lapply(seq_len(K), function(j)
          (site_offsets[j] + 1L):(site_offsets[j] + n_vec[j]/2L)))
        I2 <- unlist(lapply(seq_len(K), function(j)
          (site_offsets[j] + n_vec[j]/2L + 1L):site_offsets[j + 1L]))
        
        ## Lasso nuisance-selection on I1 EXCLUDING B columns
        nonB_idx <- setdiff(seq_len(ncol(X_noInt)), B_shift)
        X_train  <- X_noInt[I1, nonB_idx, drop = FALSE]
        Y_train  <- Y_full[I1]
        
        cv_las <- cv.glmnet(x = X_train, y = Y_train,
                            family = "gaussian", alpha = 1,
                            nfolds = 5, intercept = TRUE, standardize = TRUE)
        
        coef_all      <- as.numeric(coef(cv_las, s = "lambda.1se")) # (1 + p_nonB)
        sel_in_nb     <- which(coef_all[-1] != 0L)                  # exclude intercept
        S_hat_shift   <- if (length(sel_in_nb)) nonB_idx[sel_in_nb] else integer(0)
        
        ## Guard: ensure OLS on I2 has positive residual df
        n_I2     <- length(I2)
        max_p_ok <- max(0L, n_I2 - 1L - length(B_shift))  # intercept + |B| + |S_hat| < n_I2
        if (length(S_hat_shift) > max_p_ok) {
          betas_nb <- coef_all[-1]
          ord      <- order(abs(betas_nb[sel_in_nb]), decreasing = TRUE)
          keep     <- ord[seq_len(max_p_ok)]
          S_hat_shift <- nonB_idx[sel_in_nb[keep]]
        }
        
        ## Fit on I2: H0 uses only S_hat; H1 uses S_hat ∪ B
        X0_I2 <- cbind(1, X_noInt[I2, S_hat_shift, drop = FALSE])
        X1_I2 <- cbind(1, X_noInt[I2, union(S_hat_shift, B_shift), drop = FALSE])
        Y_I2  <- Y_full[I2]
        
        fit0_I2 <- lm.fit(X0_I2, Y_I2)
        fit1_I2 <- lm.fit(X1_I2, Y_I2)
        
        sigma0_sq_I2 <- sum(fit0_I2$residuals^2) / n_I2
        sigma1_sq_I2 <- sum(fit1_I2$residuals^2) / n_I2
        
        LR_split   <- n_I2 * log(sigma0_sq_I2 / sigma1_sq_I2)
        ## pval_split available if you want to inspect, but not stored:
        # pval_split <- 1 - pchisq(LR_split, df = length(B_shift))
        
        ## <<< CHANGED: store 0/1 rejection indicator only
        rejections_banchmark[rep] <- (LR_split > qchisq(1 - alpha, df = length(B_shift)))
        ## ── END NEW BLOCK ────────────────────────────────────────────────────
        
        ## 4·2  Null model (B deleted) ----------------------------------------
        X_H0 <- X_noInt[, -B_shift, drop = FALSE]
        
        cv0 <- cv.glmtlp(X_H0, Y_full,
                         family = "gaussian", penalty = "l0",
                         nfolds = 5, kappa = kappa_0)
        raw_supp0_shift <- which(coef(cv0)[-1] != 0)
        
        # map indices back to X_noInt positions (contiguous B block)
        first_B <- min(B_shift); len_B <- length(B_shift)
        map_back <- function(j) if (j >= first_B) j + len_B else j
        supp0    <- vapply(raw_supp0_shift, map_back, integer(1L))
        
        X0_ols   <- cbind(1, X_noInt[, supp0, drop = FALSE])
        sigma0_sq <- sum(lm.fit(X0_ols, Y_full)$residuals^2) / n_total
        
        ## 4·3  Alt model (B forced-in) ---------------------------------------
        #        pf_H1          <- rep(1, ncol(X_noInt))
        #        pf_H1[B_shift] <- 0
        
        #        cv1 <- cv.glmtlp(X_noInt, Y_full,
        #                         family = "gaussian", penalty = "l0",
        #                         nfolds = 5, penalty.factor = pf_H1,
        #                         kappa  = kappa_0)
        #        raw_supp1 <- which(coef(cv1)[-1] != 0)
        
        ## NEW: enforce nesting directly
        raw_supp1 <- sort(union(supp0, B_shift))  # \tilde B = supp0 ∪ B
        
        
        supp1     <- union(raw_supp1, B_shift)
        
        X1_ols    <- cbind(1, X_noInt[, supp1, drop = FALSE])
        sigma1_sq <- sum(lm.fit(X1_ols, Y_full)$residuals^2) / n_total
        
        ## 4·4  Likelihood-ratio test & record --------------------------------
        LRT <- n_total * log(sigma0_sq / sigma1_sq)
        z_normal   <- (LRT - len_B) / sqrt(2*len_B)
        rejections[rep]        <- (LRT > qchisq(1 - alpha, df = len_B))
        rejections_normal[rep] <- (z_normal > qnorm(1 - alpha))
      }
      
      ## 4·5  Store rejection rates for this configuration --------------------
      key                    <- paste(p, B_size, delta, sep = "_")
      results[[key]]         <- mean(rejections)
      results_normal[[key]]  <- mean(rejections_normal)
      results_banchmark[[key]] <- mean(rejections_banchmark)    # <<< CHANGED
      
      results_df_all[[row_id]] <- data.frame(
        p                       = p,
        B_size                  = B_size,
        delta                   = delta,
        rejection_rate_chisq    = results[[key]],
        rejection_rate_normal   = results_normal[[key]],
        rejection_rate_lasso_split = results_banchmark[[key]]    # <<< CHANGED
      )
      row_id <- row_id + 1L
      
      ## <<< CHANGED: removed the old "append requested metric" lines
      ## (no need to add columns after row creation anymore)
    }
  }
}

## ── 5. Collate and export ----------------------------------------------------
results_df_all <- do.call(rbind, results_df_all)
rownames(results_df_all) <- NULL
print(results_df_all)

