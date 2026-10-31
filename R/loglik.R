# Log-likelihood functions.
#
# Notation (project description, Section 2):
#   y     counts y_i
#   mu    means mu_i = exp(eta_i)
#   theta NB dispersion parameter (Var = mu + mu^2 / theta)
#   a     prior weights a_i (all 1 for an ordinary fit; 1 - tau_i in the
#         ZINB count M-step)

# NB2 log-likelihood, Section 2.1:
#   l(beta, theta) = sum_i a_i [ log Gamma(y_i + theta) - log Gamma(theta)
#                                - log y_i! + theta log theta + y_i log mu_i
#                                - (theta + y_i) log(theta + mu_i) ]
#
# The bracket is exactly log f_NB(y_i; mu_i, theta), the log of the NB2
# probability in Section 2.1, so we compute it with
# dnbinom(y, size = theta, mu = mu, log = TRUE) (allowed by Table 2).
# Why not type the formula directly? For large theta (close to the Poisson
# limit) the terms lgamma(y + theta), lgamma(theta), theta log theta and
# (theta + y) log(theta + mu) are each huge (about 1e8 when theta = 1e7) and
# nearly cancel, so the direct sum loses most of its digits to rounding.
# dnbinom() evaluates the same quantity in a cancellation-free way.
nb_loglik <- function(y, mu, theta, a = 1) {
  sum(a * stats::dnbinom(y, size = theta, mu = mu, log = TRUE))
}

# ZINB observed-data log-likelihood, Section 2.2:
#   l(beta, gamma, theta) =
#       sum_{i: y_i = 0} log[ pi_i + (1 - pi_i) f_NB(0; mu_i, theta) ]
#     + sum_{i: y_i > 0} [ log(1 - pi_i) + log f_NB(y_i; mu_i, theta) ]
# where f_NB(0; mu, theta) = (theta / (theta + mu))^theta,
# log mu_i = x_i' beta + o_i and logit pi_i = z_i' gamma.
# f_NB is computed with dnbinom() for numerical stability (see nb_loglik()).
zinb_loglik <- function(y, mu, theta, pi) {
  f0 <- stats::dnbinom(0, size = theta, mu = mu)
  ll_zero <- log(pi + (1 - pi) * f0)
  ll_pos <- log1p(-pi) + stats::dnbinom(y, size = theta, mu = mu, log = TRUE)
  sum(ifelse(y == 0, ll_zero, ll_pos))
}

# The same log-likelihood as a function of the stacked parameter vector
# par = (beta, gamma, log theta). Used for the Hessian (standard errors).
zinb_loglik_par <- function(par, X, Z, y, offset) {
  p <- ncol(X)
  q <- ncol(Z)
  beta <- par[seq_len(p)]
  gamma <- par[p + seq_len(q)]
  theta <- exp(par[p + q + 1])
  mu <- exp(drop(X %*% beta) + offset)
  pi <- stats::plogis(drop(Z %*% gamma))
  zinb_loglik(y, mu, theta, pi)
}
