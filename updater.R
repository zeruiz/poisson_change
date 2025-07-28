source("utils.R")

################# updating functions #################
update_beta <- function(beta, X, y, A, rho1, rho2, v, eta, delta, alpha) {
  grad <- grad_neg_loglik(beta, X, y)
  hess <- hess_neg_loglik(beta, X)
  
  lhs <- hess + rho1 * t(A) %*% A + rho2 * Diagonal(length(beta))
  rhs <- hess %*% beta - grad - t(A) %*% v + rho1 * t(A) %*% eta - delta + rho2 * alpha
  
  beta_new <- solve(lhs, rhs)
  return(as.numeric(beta_new))
}


update_eta <- function(beta, v, rho1, A, lambda1, gamma1) {
  a <- as.vector(A %*% beta + v / rho1)
  L <- ncol(A) / nrow(A)
  n_group <- length(a) / L
  eta <- vector("numeric", length(a))
  
  # MCP threshold for vector
  for (g in 1:n_group) {
    idx <- ((g - 1) * L + 1):(g * L)
    vec <- a[idx]
    norm_vec <- sqrt(sum(vec^2))
    if (is.na(norm_vec) || norm_vec == 0) {
      eta[idx] <- 0 # fallback safe default
    } else if (norm_vec <= gamma1 * lambda1) {
      scale <- max(0, 1 - lambda1 / (rho1 * norm_vec))
      denom <- 1 - 1 / (gamma1 * rho1)
      if (denom <= 0) { denom <- 1e-8 } # prevent divide by zero
      eta[idx] <- scale * vec / denom
    } else {
      eta[idx] <- vec
    }
  }
  
  return(matrix(eta, ncol = 1))
}

update_alpha <- function(beta, delta, rho2, lambda2, gamma2) {
  b <- beta + delta / rho2
  alpha <- numeric(length(b))
  
  # MCP threshold for scalar
  for (i in seq_along(b)) {
    abs_bi <- abs(b[i])
    denom <- 1 - 1 / (gamma2 * rho2)
    if (denom <= 0) { denom <- 1e-8 }  # prevent divide by zero
    
    if (is.na(abs_bi) || is.nan(abs_bi)) {
      alpha[i] <- 0  # fallback safe default
    } else if (abs_bi <= gamma2 * lambda2) {
      alpha[i] <- sign(b[i]) * max(0, abs_bi - lambda2 / rho2) / denom
    } else {
      alpha[i] <- b[i]
    }
  }
  return(alpha)
}

update_v <- function(v, beta, eta, A, rho1) {
  v + rho1 * (A %*% beta - eta)
}

update_delta <- function(delta, beta, alpha, rho2) {
  delta + rho2 * (beta - alpha)
}