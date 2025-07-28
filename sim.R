source("eval.R")

# concentration level i: L = 2
# replicate j: J = 2
# time points t: T = 3
# Assumption for data: i and t is complete but j may be missing
df <- data.frame('i' = c(1,1,1,1,1,1,5,5,5), 
                 'j' = c(1,1,1,2,2,2,1,1,1), 
                 't' = c(1,2,3,1,2,3,1,2,3), 
                 'y' = c(rnorm(2), rnorm(1)+3, 
                         rnorm(2)+0.5, rnorm(1)+0.5+3,
                         rnorm(2)+5, rnorm(1)+5+3))

################# test on toy example #################
L <- length(unique(df$i))
T <- length(unique(df$t))
J <- length(unique(df$j))
X <- make_X(df, L, T)
y <- sort_y(df, T)
A <- make_A(L, T)
res <- admm_optim(X, y, A, lambda1 = 0.05, lambda2 = 0.08, gamma1 = 1, gamma2 = 1)

res_pred <- get_prediction(df, res$beta, L = 2, T = 3)

estimate_replicate_difference(res_pred)
estimate_concentration_difference(res_pred)
estimate_temporal_change(res_pred)

# estimate_replicate_difference(res_pred)
# Replicate_2_minus_1 
# -0.7836702 
# 
# estimate_concentration_difference(res_pred)
# Conc_5_minus_1 
# 1.56734 
# 
# estimate_temporal_change(res_pred)
# Time_1_to_2   Time_2_to_3 
# 1.385652e+00 -2.398507e-05 