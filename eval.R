source("est.R")

################# predictions and difference check??? #################
get_prediction <- function(df, beta, L, T) {
  X <- make_X(df, L, T)
  log_lambda <- as.numeric(X %*% beta)
  
  df_sorted <- do.call(rbind, lapply(1:T, function(t) df[df$t == t, ]))
  df_sorted$log_lambda <- log_lambda
  return(df_sorted)
}

estimate_replicate_difference <- function(df_lambda) {
  agg <- aggregate(log_lambda ~ j, data = df_lambda, mean)
  diff <- diff(agg$log_lambda)
  names(diff) <- paste0("Replicate_", agg$j[2], "_minus_", agg$j[1])
  return(diff)
}

estimate_concentration_difference <- function(df_lambda) {
  agg <- aggregate(log_lambda ~ i, data = df_lambda, mean)
  diff <- diff(agg$log_lambda)
  names(diff) <- paste0("Conc_", agg$i[2], "_minus_", agg$i[1])
  return(diff)
}

estimate_temporal_change <- function(df_lambda) {
  agg <- aggregate(log_lambda ~ t, data = df_lambda, mean)
  diffs <- diff(agg$log_lambda)
  names(diffs) <- paste0("Time_", head(agg$t, -1), "_to_", tail(agg$t, -1))
  return(diffs)
}