library(ggplot2)
library(patchwork)

################# figures for the simulation study #################
# Reads results/sim_results.rds (from sim_study.R), writes figures/sim_*.pdf/png.

results <- readRDS("results/sim_results.rds")
lambda_true <- results$settings$lambda_true
L <- ncol(lambda_true)
T <- nrow(lambda_true)
dir.create("figures", showWarnings = FALSE)

method_labels <- c(BIC_1 = "BIC, C_T = 1", BIC_logT = "BIC, C_T = log T", GLM = "Unpenalized GLM")
method_colors <- c(BIC_1 = "#2a78d6", BIC_logT = "#eb6834", GLM = "#1baf7a")
scenario_labels <- c(poisson = "Poisson data", zip = "Zero-inflated data (pi0 = 0.2)")

theme_sim <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
        axis.text = element_text(color = "grey30"),
        strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "top", legend.title = element_blank())

save_fig <- function(p, name, width, height) {
  ggsave(file.path("figures", paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path("figures", paste0(name, ".png")), p, width = width, height = height, dpi = 200, bg = "white")
}

summ <- rbind(results$poisson$summary, results$zip$summary)
summ$scenario <- factor(scenario_labels[summ$scenario], levels = scenario_labels)
pen <- subset(summ, method != "GLM")

################# Figure 1: selection accuracy #################
rate_vars <- c(time_detect = "True time\nchange found", time_exact = "Time pattern\nexact",
               conc_tpr = "True conc.\nchanges found", conc_exact = "Conc. pattern\nexact")
rates <- do.call(rbind, lapply(names(rate_vars), function(v)
  aggregate(list(value = pen[[v]]), by = list(scenario = pen$scenario, method = pen$method), mean)
  |> transform(metric = rate_vars[[v]])))
rates$metric <- factor(rates$metric, levels = rate_vars)

fp_vars <- c(time_fp = "Time\n(of 5 true zeros)", conc_fp = "Concentration\n(of 31 true zeros)")
fps <- do.call(rbind, lapply(names(fp_vars), function(v)
  aggregate(list(value = pen[[v]]), by = list(scenario = pen$scenario, method = pen$method), mean)
  |> transform(metric = fp_vars[[v]])))
fps$metric <- factor(fps$metric, levels = fp_vars)

bars <- function(d, ylab, ylim = NULL) {
  ggplot(d, aes(metric, value, fill = method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7, color = "white", linewidth = 0.5) +
    geom_text(aes(label = formatC(value, format = "f", digits = 2)), position = position_dodge(width = 0.75),
              vjust = -0.4, size = 2.8, color = "grey20") +
    facet_wrap(~ scenario) +
    scale_fill_manual(values = method_colors, labels = method_labels) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = ylab) + theme_sim
}
p1 <- bars(rates, "Proportion of datasets", c(0, 1.08)) +
  bars(fps, "Mean number of false positives") + theme(legend.position = "none") +
  plot_layout(ncol = 1, heights = c(1.2, 1))
save_fig(p1, "sim_detection", 9, 7.5)

################# Figure 2: estimated mean counts over days #################
traj <- do.call(rbind, lapply(names(scenario_labels), function(sc) {
  do.call(rbind, lapply(c("BIC_logT", "GLM"), function(m) {
    arr <- simplify2array(lapply(results[[sc]]$mu_hat, function(x) exp(x[[m]]))) # T x L x rep
    grid <- expand.grid(t = 1:T, i = 1:L)
    grid$median <- apply(arr, c(1, 2), median)[cbind(grid$t, grid$i)]
    grid$lo <- apply(arr, c(1, 2), quantile, 0.1)[cbind(grid$t, grid$i)]
    grid$hi <- apply(arr, c(1, 2), quantile, 0.9)[cbind(grid$t, grid$i)]
    transform(grid, method = m, scenario = scenario_labels[[sc]],
              truth = ((1 - results$settings$scenarios[[sc]]) * lambda_true)[cbind(t, i)])
  }))
}))
traj$scenario <- factor(traj$scenario, levels = scenario_labels)
traj$level <- factor(paste("Level", traj$i))
dodge <- position_dodge(width = 0.35)

p2 <- ggplot(traj, aes(t, median, color = method)) +
  geom_step(aes(y = truth, linetype = "Truth"), color = "grey20", direction = "mid", linewidth = 0.5) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.6, position = dodge) +
  geom_point(size = 1.8, position = dodge) +
  facet_grid(scenario ~ level) +
  scale_color_manual(values = method_colors, labels = method_labels) +
  scale_linetype_manual(values = c(Truth = "dashed")) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
  scale_x_continuous(breaks = 1:T) +
  labs(x = "Day", y = "Estimated mean count (log scale)",
       caption = "Points: median over replicates; bars: 10th-90th percentile.") +
  theme_sim
save_fig(p2, "sim_trajectories", 11, 5.5)

################# Figure 3: estimation error #################
p3 <- ggplot(summ, aes(method, rmse_mean, fill = method)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8, color = "grey25", linewidth = 0.4) +
  facet_wrap(~ scenario) +
  scale_fill_manual(values = method_colors, labels = method_labels) +
  scale_x_discrete(labels = method_labels) +
  labs(x = NULL, y = "RMSE of estimated mean counts") +
  theme_sim + theme(legend.position = "none")
save_fig(p3, "sim_rmse", 8, 4)
