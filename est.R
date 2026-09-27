source("updater.R")

################# ADMM optimizer #################
# MCP proximal updates need gamma * rho > 1 (see Eq. sol_eta / sol_alpha).
# Stops when the primal residuals (A beta - eta, beta - alpha) and the change
# in beta are all below tol.
admm_optim <- function(X, y, A, lambda1, lambda2, gamma1, gamma2,
                       rho1 = 1, rho2 = 1, max_iter = 500, tol = 1e-4) {
  if (gamma1 * rho1 <= 1) stop("need gamma1 * rho1 > 1 for the MCP update of eta")
  if (gamma2 * rho2 <= 1) stop("need gamma2 * rho2 > 1 for the MCP update of alpha")

  beta <- rep(0, ncol(X))
  eta <- A %*% beta
  alpha <- beta
  v <- rep(0, length(eta))
  delta <- rep(0, length(beta))
  obj_Q <- rep(NA_real_, max_iter + 1)
  obj_Q[1] <- get_Q(beta, X, y, A, lambda1, lambda2)
  converged <- FALSE

  for (iter in 1:max_iter) {
    beta_new <- update_beta(beta, X, y, A, rho1, rho2, v, eta, delta, alpha)
    eta <- update_eta(beta_new, v, rho1, A, lambda1, gamma1)
    alpha <- update_alpha(beta_new, delta, rho2, lambda2, gamma2)
    v <- update_v(v, beta_new, eta, A, rho1)
    delta <- update_delta(delta, beta_new, alpha, rho2)

    obj_Q[iter + 1] <- get_Q(beta_new, X, y, A, lambda1, lambda2)

    r_eta <- sqrt(sum(as.numeric(A %*% beta_new - eta)^2))
    r_alpha <- sqrt(sum((beta_new - alpha)^2))
    d_beta <- sqrt(sum((beta_new - beta)^2))
    beta <- beta_new
    if (max(r_eta, r_alpha, d_beta) < tol) {
      converged <- TRUE
      break
    }
  }

  list(beta = beta, eta = eta, alpha = alpha, v = v, delta = delta,
       obj_Q = obj_Q[1:(iter + 1)], iter = iter, converged = converged)
}
