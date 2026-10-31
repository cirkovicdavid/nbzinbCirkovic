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
