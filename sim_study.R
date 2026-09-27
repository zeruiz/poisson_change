source("eval.R")
library(parallel)

################# simulation study #################
# Setting (main.tex, Section "Simulation"): L = 6 concentration levels, T = 7 days,
# n_it drawn from {8, 9, 10}. Mean counts
#   lambda_t = (0.1, ..., 0.1)          for t = 1, 2, 3
#   lambda_t = (10, 10, 10, 10, 5, 5)   for t = 4, ..., 7
# i.e. mu_t = log(lambda_t): one time change (day 3 -> 4) and, from day 4 on,
# one concentration change (level 4 -> 5).
# Scenarios: "poisson" (the model is correct) and "zip" (zero-inflated Poisson
# with structural-zero probability pi0; the fitted model is still Poisson, so the
# target is the marginal mean (1 - pi0) * lambda, which has the same change pattern).
#
# Usage: Rscript sim_study.R [n_rep]   (default 100)

L <- 6
T <- 7
lambda_true <- rbind(matrix(0.1, 3, L),
                     matrix(rep(c(10, 10, 10, 10, 5, 5), each = 4), 4, L))
scenarios <- list(poisson = 0, zip = 0.2) # pi0

true_time <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE) # change between day t and t + 1
true_conc <- matrix(FALSE, T, L)                        # beta_it != 0, levels i >= 2
true_conc[4:7, 5] <- TRUE

grid <- expand.grid(lambda1 = c(0.01, 0.02, 0.05, 0.1, 0.2),
                    lambda2 = c(0.01, 0.02, 0.05, 0.1, 0.2))
gamma <- 3
mu_min <- log(0.01) # lower bound on mu (see est.R); to be discussed
max_iter <- 2000
C_T <- c(BIC_1 = 1, BIC_logT = log(T))

gen_data <- function(lambda, pi0) {
  do.call(rbind, lapply(1:T, function(t) do.call(rbind, lapply(1:L, function(i) {
    n <- sample(8:10, 1)
    y <- rpois(n, lambda[t, i]) * rbinom(n, 1, 1 - pi0)
    data.frame(i = i, j = 1:n, t = t, y = y)
  }))))
}

# Distinct nonzero parameters: time segments are split where eta_t != 0, and within
# a segment count the levels whose alpha is nonzero at any of its days.
count_df <- function(eta_nz, alpha_nz) {
  seg <- cumsum(c(1, eta_nz))
  sum(sapply(unique(seg), function(s) sum(colSums(alpha_nz[seg == s, , drop = FALSE]) > 0)))
}

fit_one <- function(X, y, A, lambda1, lambda2) {
  res <- admm_optim(X, y, A, lambda1, lambda2, gamma, gamma,
                    max_iter = max_iter, mu_min = mu_min)
  H <- kronecker(Diagonal(T), make_H(L))
  eta_norm <- sapply(1:(T - 1), function(t) sqrt(sum(res$eta[((t - 1) * L + 1):(t * L)]^2)))
  alpha_nz <- matrix(res$alpha != 0, T, L, byrow = TRUE)
  list(mu_hat = matrix(as.numeric(H %*% res$beta), T, L, byrow = TRUE),
       eta_nz = eta_norm > 0, alpha_nz = alpha_nz,
       nll = neg_loglik(res$beta, X, y),
       df = count_df(eta_norm > 0, alpha_nz),
       converged = res$converged, iter = res$iter)
}

evaluate <- function(fit, target) {
  a <- fit$alpha_nz[, -1]
  tc <- true_conc[, -1]
  c(time_detect = fit$eta_nz[3],
    time_fp = sum(fit$eta_nz[-3]),
    time_exact = all(fit$eta_nz == true_time),
    conc_tpr = mean(a[tc]),
    conc_fp = sum(a & !tc),
    conc_exact = all(a == tc),
    rmse_log = sqrt(mean((fit$mu_hat - log(target))^2)),
    rmse_mean = sqrt(mean((exp(fit$mu_hat) - target)^2)))
}

run_rep <- function(rep, pi0) {
  set.seed(1000 + rep)
  df <- gen_data(lambda_true, pi0)
  target <- (1 - pi0) * lambda_true
  X <- make_X(df, L, T)
  y <- sort_y(df, T)
  A <- make_A(L, T)
  n <- length(y)

  fits <- lapply(seq_len(nrow(grid)), function(k) fit_one(X, y, A, grid$lambda1[k], grid$lambda2[k]))
  nll <- sapply(fits, `[[`, "nll")
  dfs <- sapply(fits, `[[`, "df")

  out <- lapply(names(C_T), function(crit) {
    k <- which.min(2 * nll + C_T[[crit]] * dfs * log(n) / n)
    f <- fits[[k]]
    list(summary = data.frame(rep = rep, method = crit, t(evaluate(f, target)),
                              lambda1 = grid$lambda1[k], lambda2 = grid$lambda2[k],
                              df = f$df, converged = f$converged, iter = f$iter),
         mu_hat = f$mu_hat)
  })

  # unpenalized Poisson GLM with one mean per cell = log of the cell average
  ybar <- tapply(df$y, list(df$t, df$i), mean)
  mu_glm <- pmax(log(ybar), mu_min)
  glm <- data.frame(rep = rep, method = "GLM",
                    rmse_log = sqrt(mean((mu_glm - log(target))^2)),
                    rmse_mean = sqrt(mean((exp(mu_glm) - target)^2)))

  list(summary = dplyr::bind_rows(lapply(out, `[[`, "summary"), glm),
       mu_hat = c(setNames(lapply(out, `[[`, "mu_hat"), names(C_T)), list(GLM = mu_glm)),
       conv_all = mean(sapply(fits, `[[`, "converged")))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  n_rep <- if (length(args) > 0) as.integer(args[1]) else 100
  dir.create("results", showWarnings = FALSE)

  results <- lapply(names(scenarios), function(sc) {
    t0 <- Sys.time()
    reps <- mclapply(1:n_rep, run_rep, pi0 = scenarios[[sc]],
                     mc.cores = max(1, detectCores() - 1))
    cat(sprintf("%s: %d reps in %.1f min\n", sc, n_rep,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    list(summary = cbind(scenario = sc, dplyr::bind_rows(lapply(reps, `[[`, "summary"))),
         mu_hat = lapply(reps, `[[`, "mu_hat"),
         conv_all = sapply(reps, `[[`, "conv_all"))
  })
  names(results) <- names(scenarios)
  results$settings <- list(lambda_true = lambda_true, scenarios = scenarios, grid = grid,
                           gamma = gamma, mu_min = mu_min, max_iter = max_iter, C_T = C_T)
  saveRDS(results, "results/sim_results.rds")

  summ <- dplyr::bind_rows(results$poisson$summary, results$zip$summary)
  tab <- dplyr::summarise(dplyr::group_by(summ, scenario, method),
                          dplyr::across(c(time_detect:rmse_mean, converged), ~ mean(.x)),
                          .groups = "drop")
  write.csv(tab, "results/sim_summary.csv", row.names = FALSE)
  print(as.data.frame(tab), digits = 3)
  for (sc in names(scenarios))
    cat(sprintf("%s: share of all grid fits converged = %.3f\n", sc, mean(results[[sc]]$conv_all)))
}
