source("updater.R")

################# ADMM optimizer #################
admm_optim <- function(X, y, A, lambda1, lambda2, gamma1, gamma2,
                       rho1 = 1, rho2 = 1, max_iter = 100, tol = 1e-4) {
  beta <- rep(0, ncol(X))
  eta <- A %*% beta
  alpha <- beta
  v <- rep(0, length(eta))
  delta <- rep(0, length(beta))
  obj_Q <- rep(0, max_iter + 1)
  obj_Q[1] <- get_Q(beta, X, y, A, lambda1, lambda2)
  
  for (iter in 1:max_iter) {
    beta_new <- update_beta(beta, X, y, A, rho1, rho2, v, eta, delta, alpha)
    eta_new <- update_eta(beta_new, v, rho1, A, lambda1, gamma1)
    alpha_new <- update_alpha(beta_new, delta, rho2, lambda2, gamma2)
    v <- update_v(v, beta_new, eta_new, A, rho1)
    delta <- update_delta(delta, beta_new, alpha_new, rho2)
    
    obj_Q[iter + 1] <- get_Q(beta_new, X, y, A, lambda1, lambda2)
    
    if (sqrt(sum((beta_new - beta)^2)) < tol) {break}
    beta <- beta_new
    eta <- eta_new
    alpha <- alpha_new
  }
  
  list(beta = beta, eta = eta, alpha = alpha, v = v, delta = delta, obj_Q = obj_Q)
}