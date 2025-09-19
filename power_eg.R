## -- 0.  Local library path ----------------------------------------
.libPaths(c("/panfs/jay/groups/27/shenx/zhao1118/distributed_algorithm/R_libs",
            .libPaths()))
library(MASS)      # recommended, usually pre-installed
library(glmtlp)    # installed in your private lib
num_parallel  <- 20          # SLURM tasks (array length)


## -- 1.  Global grids & replication bookkeeping --------------------
p_grid        <- c(50, 100, 500, 1000)
B_sizes       <- c(1, 5, 10)
delta_vec     <- c(0.0, 0.1, 0.2, 0.3)

n_rep_total   <- 1000        # desired reps / parameter triple
reps_per_job  <- n_rep_total / num_parallel   # 50

## -- 1·1  chunk id from SLURM --------------------------------------
rep_chunk <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
if (is.na(rep_chunk) || rep_chunk < 1 || rep_chunk > num_parallel)
  stop("SLURM_ARRAY_TASK_ID must be 1–", num_parallel, call. = FALSE)

## deterministic seed stream: unique base per chunk
base_seed <- 74000L + rep_chunk * 100000L
set.seed(base_seed)

## -- 2.  Fixed simulation settings ---------------------------------
K            <- 3
n_vec        <- rep(100, K)
p_vec        <- c(10, 10, 10)
rho_eq       <- 0.1
sigma_error  <- 1

length_b0    <- 3                       # global slopes (1,2,3)
true_sparse  <- length_b0 + 3 - 1
start_k      <- ceiling((true_sparse - 1) / 2)
end_k        <- ceiling(true_sparse * 1.20)
kappa_0      <- sort(rep(seq.int(start_k, end_k), length.out = 60))

## -- 3.  Helper functions  (unchanged from earlier) -----------------

generate_design <- function(K, n, p0, pj, rho) {
  max_p <- max(p0, pj)
  Sigma <- matrix(rho, max_p, max_p); diag(Sigma) <- 1
  lapply(seq_len(K), function(j) {
    X <- if (p0 > 1) {
      cbind(Intercept = 1,
            MASS::mvrnorm(n, rep(0, p0 - 1), Sigma[1:(p0 - 1), 1:(p0 - 1)]))
    } else matrix(1, n, 1, dimnames = list(NULL, "Intercept"))
    W <- MASS::mvrnorm(n, rep(0, pj), Sigma[1:pj, 1:pj])
    list(X = X, W = W)
  })
}

generate_beta0 <- function(p0, B_idx, delta, length_b0= 3) {
  b0 <- numeric(p0); b0[1:3] <- c(1, 2, 3); b0[B_idx[1]] <- delta; b0
}
generate_beta_list <- function(p_vec) {
  lapply(seq_along(p_vec), \(j) { b <- numeric(p_vec[j]); b[1] <- j; b })
}

generate_data <- function(design, beta0, beta_list, sigma = 1) {
  mapply(function(site, bj) {
    eps <- rnorm(nrow(site$X), sd = sigma)
    drop(site$X %*% beta0 + site$W %*% bj + eps)
  }, design, beta_list, SIMPLIFY = FALSE)
}

assemble_full_model <- function(design, beta0, beta_list, Y_list) {
  W_list <- lapply(design, `[[`, "W")
  p_vec  <- sapply(W_list, ncol)
  rows <- lapply(seq_along(design), function(j) {
    Xj <- design[[j]]$X; Wj <- W_list[[j]]; nj <- nrow(Xj)
    left  <- if (j > 1) matrix(0, nj, sum(p_vec[1:(j - 1)])) else NULL
    right <- if (j < length(design)) matrix(0, nj, sum(p_vec[(j + 1):length(design)])) else NULL
    cbind(Xj, left, Wj, right)
  })
  list(
    design   = do.call(rbind, rows),
    beta     = c(beta0, unlist(beta_list)),
    response = unlist(Y_list)
  )
}


## -- 4.  Output directory ------------------------------------------
out_dir <- "/panfs/jay/groups/27/shenx/zhao1118/distributed_algorithm/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## -- 5.  Main three-level loop -------------------------------------
combo_counter <- 0L   # for seeding uniqueness

for (p in p_grid) {
  p0 <- p - sum(p_vec)                # global slopes excl. intercept
  stopifnot(p0 > length_b0)           # sanity
  
  for (B_size in B_sizes) {
    B_idx   <- (length_b0 + 1):(length_b0 + B_size)
    B_shift <- B_idx - 1              # index inside X_noInt
    len_B   <- B_size
    
    for (delta in delta_vec) {
      
      combo_counter <- combo_counter + 1L
      combo_seed    <- base_seed + combo_counter * 10^6L  # huge stride
      
      ## ---- 5·1  replication loop for this triple -----------------
      for (r in seq_len(reps_per_job)) {
        
        set.seed(combo_seed + r)      # unique & reproducible
        
        ## 5·1·1  Data generation
        design_list <- lapply(seq_len(K), function(j)
          generate_design(1, n_vec[j], p0, p_vec[j], rho_eq)[[1]]
        )
        beta0_true  <- generate_beta0(p0, B_idx, delta, length_b0)
        beta_list   <- generate_beta_list(p_vec)
        Y_list      <- generate_data(design_list, beta0_true,
                                     beta_list, sigma_error)
        
        full        <- assemble_full_model(design_list, beta0_true,
                                           beta_list, Y_list)
        X_full      <- full$design
        Y_full      <- full$response
        n_total     <- length(Y_full)
        X_noInt     <- X_full[, -1, drop = FALSE]
        
        ## 5·1·2  H0 (B removed)
        X_H0 <- X_noInt[, -B_shift, drop = FALSE]
        cv0  <- cv.glmtlp(X_H0, Y_full, family = "gaussian",
                          penalty = "l0", nfolds = 5, kappa = kappa_0)
        
        raw_supp0_shift <- which(coef(cv0)[-1] != 0)
        first_B <- min(B_shift)
        supp0   <- if (length(raw_supp0_shift))
          raw_supp0_shift + (raw_supp0_shift >= first_B) * len_B
        else integer(0)
        X0_ols   <- cbind(1, X_noInt[, supp0, drop = FALSE])
        sigma0_sq<- sum(lm.fit(X0_ols, Y_full)$res^2) / n_total
        
        ## 5·1·3  H1 (B forced-in)
        pf_H1 <- rep(1, ncol(X_noInt)); pf_H1[B_shift] <- 0
        cv1 <- cv.glmtlp(X_noInt, Y_full, family = "gaussian",
                         penalty = "l0", nfolds = 5,
                         penalty.factor = pf_H1, kappa = kappa_0)
        
        raw_supp1 <- which(coef(cv1)[-1] != 0)
        supp1     <- union(raw_supp1, B_shift)
        X1_ols    <- cbind(1, X_noInt[, supp1, drop = FALSE])
        sigma1_sq <- sum(lm.fit(X1_ols, Y_full)$res^2) / n_total
        
        ## 5·1·4  LRT
        LRT      <- n_total * log(sigma0_sq / sigma1_sq)
        decision <- (LRT > qchisq(0.95, df = len_B))
        
        ## 5·1·5  Save one replicate
        result_i <- list(
          p         = p,
          B_size    = B_size,
          delta     = delta,
          rep_chunk = rep_chunk,
          repl_id   = r,
          reject    = decision,
          metric    = if (delta == 0) "size" else "power"
        )
        
        file_tag_i <- sprintf("p%d_B%d_delta%.1f_rep%02d_iter%03d",
                              p, B_size, delta, rep_chunk, r)
        
        saveRDS(result_i,
                file = file.path(out_dir, paste0(file_tag_i, ".RDS")))
      } # end replication loop
    }   # end delta loop
  }     # end B_size loop
}       # end p grid loop
