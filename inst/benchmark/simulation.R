# Simulation study (project description, Section 5.1): does zinbreg()
# recover known parameters?
#
# Design: R = 500 data sets of size n = 500 from a ZINB model with one
# covariate in each part,
#   count:  log mu_i  = beta0 + beta1 x_i,     x_i ~ N(0, 1)
#   zero:   logit pi_i = gamma0 + gamma1 z_i,  z_i ~ N(0, 1)
#   theta = 2
# For each parameter we report the bias, the mean of the estimated standard
# errors vs. the Monte Carlo standard deviation of the estimates, and the
# empirical coverage of 95% Wald intervals. Results are saved to
# inst/benchmark/results/simulation.rds and used by the report.
#
# Run from the package root:  Rscript inst/benchmark/simulation.R

library(nbzinbCirkovic)

set.seed(20261031)
R <- 500
n <- 500
truth <- c(beta0 = 1, beta1 = 0.5, gamma0 = -1, gamma1 = 0.8,
           log_theta = log(2))

simulate_zinb <- function(n, truth) {
  x <- rnorm(n)
  z <- rnorm(n)
  mu <- exp(truth[["beta0"]] + truth[["beta1"]] * x)
  pi <- plogis(truth[["gamma0"]] + truth[["gamma1"]] * z)
  structural <- rbinom(n, 1, pi)
  y <- ifelse(structural == 1, 0,
              rnbinom(n, size = exp(truth[["log_theta"]]), mu = mu))
  data.frame(y = y, x = x, z = z)
}

est <- se <- matrix(NA_real_, R, length(truth),
                    dimnames = list(NULL, names(truth)))
converged <- logical(R)
prop_zero <- numeric(R)
elapsed <- system.time(
  for (r in seq_len(R)) {
    d <- simulate_zinb(n, truth)
    prop_zero[r] <- mean(d$y == 0)
    fit <- suppressWarnings(zinbreg(y ~ x | z, data = d))
    converged[r] <- fit$converged
    est[r, ] <- c(coef(fit), log(fit$theta))
    se[r, ] <- c(sqrt(diag(vcov(fit))), fit$SE.logtheta)
  }
)[["elapsed"]]

lower <- est - qnorm(0.975) * se
upper <- est + qnorm(0.975) * se
covered <- sweep(lower, 2, truth, "<=") & sweep(upper, 2, truth, ">=")

summary_tab <- data.frame(
  parameter = names(truth),
  true = truth,
  mean_est = colMeans(est),
  bias = colMeans(est) - truth,
  mean_se = colMeans(se),
  mc_sd = apply(est, 2, sd),
  coverage = colMeans(covered),
  row.names = NULL
)
print(summary_tab, digits = 3)
cat("Converged:", sum(converged), "of", R, " Mean proportion of zeros:",
    round(mean(prop_zero), 3), " Time (s):", round(elapsed, 1), "\n")

saveRDS(list(summary = summary_tab, est = est, se = se, truth = truth,
             R = R, n = n, converged = converged, prop_zero = prop_zero,
             elapsed = elapsed, seed = 20261031),
        file.path("inst", "benchmark", "results", "simulation.rds"))
