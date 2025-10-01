#!/usr/bin/env Rscript
.libPaths()

library(MASS)
library(glmtlp)
library(hdi)

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
  b0[1:length_b0] <- seq_len(length_b0)
  b0[B_idx[1]]   <- delta
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
p_grid       <- c(50)
B_sizes      <- c(1, 5, 10)
delta_vec    <- c(0.0, 0.1, 0.2, 0.3)
n_vec        <- rep(100, 3)
p_vec        <- c(10, 10, 10)
rho_eq       <- 0.5
n_rep        <- 1000
sigma_error  <- 1

alpha        <- 0.05
length_b0    <- 3
true_sparse  <- length_b0 + 3 - 1
kappa_0      <- as.integer(rep(true_sparse, 30))

K  <- length(n_vec)

## ── 3. Containers for results ───────────────────────────────────────────────
results         <- list()
results_normal  <- list()
results_de_lass <- list()
results_df_all  <- list()
row_id          <- 1L

timing <- list(cv0 = 0, cv1 = 0, lasso_proj = 0)

time_start <- Sys.time()

max_runs_env <- Sys.getenv("MAX_RUNS", NA_character_)
max_runs <- suppressWarnings(as.integer(max_runs_env))
if (is.na(max_runs) || max_runs < 1L) max_runs <- Inf
run_counter <- 0L
early_stop  <- FALSE

## ── 4. Main simulation loops ────────────────────────────────────────────────
for (p in p_grid) {
  p0 <- p - sum(p_vec)

  for (B_size in B_sizes) {
    B_idx   <- (length_b0 + 1):(length_b0 + B_size)
    B_shift <- B_idx - 1

    for (delta in delta_vec) {

      rejections          <- logical(n_rep)
      rejections_normal   <- logical(n_rep)
      pvals_de_lasso      <- numeric(n_rep)

      for (rep in seq_len(n_rep)) {

        ## 4·1  Data generation ----------------------------------------------
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

        ## 4·2  Null model ---------------------------------------------------
        X_H0 <- X_noInt[, -B_shift, drop = FALSE]

        t_cv0 <- proc.time()
        cv0 <- cv.glmtlp(X_H0, Y_full,
                         family = "gaussian", penalty = "l0",
                         nfolds = 5, kappa = kappa_0)
        timing$cv0 <- timing$cv0 + (proc.time() - t_cv0)["elapsed"]
        raw_supp0_shift <- which(coef(cv0)[-1] != 0)

        first_B <- min(B_shift); len_B <- length(B_shift)
        map_back <- function(j) if (j >= first_B) j + len_B else j
        supp0    <- vapply(raw_supp0_shift, map_back, integer(1L))

        X0_ols   <- cbind(1, X_noInt[, supp0, drop = FALSE])
        sigma0_sq <- sum(lm.fit(X0_ols, Y_full)$residuals^2) / n_total

        ## 4·3  Alt model ----------------------------------------------------
        pf_H1          <- rep(1, ncol(X_noInt))
        pf_H1[B_shift] <- 0

        t_cv1 <- proc.time()
        cv1 <- cv.glmtlp(X_noInt, Y_full,
                         family = "gaussian", penalty = "l0",
                         nfolds = 5, penalty.factor = pf_H1,
                         kappa  = kappa_0)
        timing$cv1 <- timing$cv1 + (proc.time() - t_cv1)["elapsed"]
        raw_supp1 <- which(coef(cv1)[-1] != 0)
        supp1     <- union(raw_supp1, B_shift)

        X1_ols    <- cbind(1, X_noInt[, supp1, drop = FALSE])
        sigma1_sq <- sum(lm.fit(X1_ols, Y_full)$residuals^2) / n_total

        ## 4·4  Likelihood-ratio tests --------------------------------------
        LRT <- n_total * log(sigma0_sq / sigma1_sq)
        z_normal   <- (LRT - len_B) / sqrt(2 * len_B)
        rejections[rep]        <- (LRT > qchisq(0.95, df = len_B))
        rejections_normal[rep] <- (z_normal > qnorm(0.95))

        ## 4·5  Debiased Lasso group p-value --------------------------------
        t_lasso_proj <- proc.time()
        fit_dl <- hdi::lasso.proj(
          x = X_noInt,
          y = Y_full,
          standardize = TRUE,
          suppress.grouptesting = FALSE
        )
        timing$lasso_proj <- timing$lasso_proj + (proc.time() - t_lasso_proj)["elapsed"]

        p_group <- tryCatch(
          as.numeric(fit_dl$groupTest(B_shift)),
          error = function(e) NA_real_
        )
        if (is.na(p_group) || length(p_group) != 1L) p_group <- 1
        pvals_de_lasso[rep] <- p_group
      }

      ## 4·6  Store results -------------------------------------------------
      key                         <- paste(p, B_size, delta, sep = "_")
      results[[key]]              <- mean(rejections)
      results_normal[[key]]       <- mean(rejections_normal)
      results_de_lass[[key]]      <- mean(pvals_de_lasso < alpha)

      results_df_all[[row_id]] <- data.frame(
        p                        = p,
        B_size                   = B_size,
        delta                    = delta,
        rejection_rate_chisq     = results[[key]],
        rejection_rate_normal    = results_normal[[key]],
        rejection_rate_deblasso  = results_de_lass[[key]],
        mean_pvalue_deblasso     = mean(pvals_de_lasso)
      )
      row_id <- row_id + 1L

      run_counter <- run_counter + 1L
      if (is.finite(max_runs) && run_counter >= max_runs) {
        early_stop <- TRUE
        break
      }
    }
    if (early_stop) break
  }
  if (early_stop) break
}

## ── 5. Collate and export ----------------------------------------------------
results_df_all <- do.call(rbind, results_df_all)
rownames(results_df_all) <- NULL
print(results_df_all)

saveRDS(results_df_all, "simulation_results_deblasso.rds")

time_end <- Sys.time()
elapsed_secs <- as.numeric(difftime(time_end, time_start, units = "secs"))
cat(sprintf("Total elapsed time: %.2f seconds\n", elapsed_secs))
cat(sprintf("Cumulative timing (seconds) — cv.glmtlp H0: %.2f, cv.glmtlp H1: %.2f, lasso.proj: %.2f\n",
            timing$cv0, timing$cv1, timing$lasso_proj))

