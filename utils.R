library(Matrix)

################# utils #################
make_H <- function(L) {
  H <- matrix(0, nrow = L, ncol = L)
  H[lower.tri(H, diag = TRUE)] <- 1
  return(H)
}

make_X <- function(df, L, T) {
  df$i_mapped <- as.numeric(factor(df$i)) # when the difference between concentration level is not 1 
  H <- make_H(L)
  X_blocks <- vector("list", T)
  for (t in 1:T) {
    df_t <- df[df$t == t, ]
    R_t <- model.matrix(~ factor(i_mapped, levels = 1:L) - 1, data = df_t)
    R_t <- Matrix(R_t, sparse = TRUE)
    X_blocks[[t]] <- R_t %*% H
  }
  X <- bdiag(X_blocks)
  return(X)
}

make_A <- function(L, T) {
  M <- bandSparse(T - 1, T, k = 0:1, diagonals = list(rep(-1, T - 1), rep(1, T - 1)))
  A <- kronecker(M, Diagonal(L))
  return(A)
}

# This function is used when df's layout is i, j, t
# but we actually needs y_it in blocks
sort_y <- function(df, T) {
  y_list <- vector("list", T)
  for (t in 1:T) {
    df_t <- df[df$t == t, ]
    y_list[[t]] <- df_t$y
  }
  y <- unlist(y_list)
  return(y)
}

################# objective function #################
neg_loglik <- function(beta, X, y) {
  log_lambda <- X %*% beta
  mu <- exp(log_lambda)
  -mean(y * log_lambda - mu)
}

grad_neg_loglik <- function(beta, X, y) {
  log_lambda <- X %*% beta
  mu <- exp(log_lambda)
  -t(X) %*% (y - mu) / length(y)
}

hess_neg_loglik <- function(beta, X) {
  log_lambda <- X %*% beta
  mu <- exp(log_lambda)
  V <- Diagonal(x = as.numeric(mu))
  t(X) %*% V %*% X / nrow(X)
}

get_Q <- function(beta, X, y, A, lambda1, lambda2) {
  nll <- neg_loglik(beta, X, y)
  
  L <- ncol(A) / nrow(A)
  Abeta <- A %*% beta
  group_norms <- sapply(seq(1, length(Abeta), by = L), function(k) {
    sqrt(sum(Abeta[k:(k + L - 1)]^2))
  })
  group_penalty <- lambda1 * sum(group_norms)
  l1_penalty <- lambda2 * sum(abs(beta))
  
  nll + group_penalty + l1_penalty
}