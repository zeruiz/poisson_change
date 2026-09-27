# poisson_change

Change detection over time and concentration level for Poisson counts:
log(lambda_it) = mu_it, with a group MCP penalty on beta_t - beta_{t-1} and an MCP
penalty on beta_it (differences across concentration levels), fit by ADMM.

| File | Contents |
|---|---|
| `utils.R` | design matrices, likelihood, objective, `clamp_mu` |
| `updater.R` | ADMM updates for beta, eta, alpha and the dual variables |
| `est.R` | `admm_optim` |
| `eval.R` | predictions and summary differences |
| `sim.R` | toy example |
| `sim_study.R` | simulation study (writes `results/`) |
| `plot_sim.R` | simulation figures (writes `figures/`) |

Reproduce the simulation study (about 35 minutes on 10 cores):

```
Rscript sim_study.R 100
Rscript plot_sim.R
```

Open issue: when every replicate in a cell is zero, the likelihood has no finite
maximizer and MCP does not bound it, so the simulation uses a lower bound
`mu_min = log(0.01)` (`admm_optim(..., mu_min = )`). To be discussed.
