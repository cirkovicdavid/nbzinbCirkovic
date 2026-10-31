# Inner algorithms for the NB model (project description, Section 3.1).
# These are internal (not exported). nb_fit() in nbreg.R alternates them, and
# the ZINB count M-step will reuse them with prior weights a_i = 1 - tau_i.

# Weighted least squares: returns the beta that solves the normal equations
#   (X' W X) beta = X' W z,   W = diag(w).
# Section 3.1(a) allows solve(), qr() or chol(). We use the QR decomposition
# of W^(1/2) X because it is numerically more stable than forming X'WX:
# minimizing ||W^(1/2) z - W^(1/2) X beta||^2 gives the same beta.
# (X * sw multiplies row i of X by sw[i].)
wls_solve <- function(X, z, w) {
  sw <- sqrt(w)
  qrx <- qr(X * sw)
  if (qrx$rank < ncol(X)) {
    stop("the model matrix is rank deficient; ",
         "remove redundant (collinear) covariates", call. = FALSE)
  }
  beta <- qr.coef(qrx, z * sw)
  names(beta) <- colnames(X)
  beta
}

# IRLS working weights, Section 3.1(a):
#   w_i = a_i (d mu_i / d eta_i)^2 / V(mu_i) = a_i theta mu_i / (theta + mu_i),
# because d mu / d eta = mu for the log link and V(mu) = mu + mu^2 / theta.
# theta = Inf is the Poisson limit, V(mu) = mu, which gives w_i = a_i mu_i.
# We use that case once, for the Poisson starting fit (Section 3.1(c)).
irls_weights <- function(mu, theta, a) {
  if (is.infinite(theta)) a * mu else a * theta * mu / (theta + mu)
}

# IRLS for beta with theta held fixed, Section 3.1(a).
#   X      n x p model matrix
#   y      counts
#   offset o_i (zeros if there is no offset)
#   a      prior weights
#   theta  fixed dispersion (Inf = Poisson)
#   beta   starting value; NULL means start from mu_i = y_i + 0.1
#          (the +0.1 avoids log(0) when y_i = 0)
#   tol    stop when max |beta(t+1) - beta(t)| < tol
irls_nb <- function(X, y, offset, a, theta, beta = NULL, tol = 1e-8,
                    maxit = 100) {
  if (is.null(beta)) {
    mu <- y + 0.1
    eta <- log(mu)
  } else {
    eta <- drop(X %*% beta) + offset
    mu <- exp(eta)
  }
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    # Working weights and working response, Section 3.1(a):
    #   w_i = a_i theta mu_i / (theta + mu_i)
    #   z_i = (eta_i - o_i) + (y_i - mu_i) / mu_i
    w <- irls_weights(mu, theta, a)
    z <- (eta - offset) + (y - mu) / mu
    # beta(t+1) = (X' W X)^(-1) X' W z
    beta_new <- wls_solve(X, z, w)
    eta <- drop(X %*% beta_new) + offset
    mu <- exp(eta)
    if (any(!is.finite(mu))) {
      stop("IRLS diverged (non-finite fitted means)", call. = FALSE)
    }
    if (!is.null(beta) && max(abs(beta_new - beta)) < tol) {
      beta <- beta_new
      converged <- TRUE
      break
    }
    beta <- beta_new
  }
  list(beta = beta, eta = eta, mu = mu, iter = iter, converged = converged)
}

# Score and second derivative of the NB log-likelihood in theta, with mu
# fixed, Section 3.1(b):
#   U(theta) = sum_i a_i [ psi(y_i + theta) - psi(theta) + log theta + 1
#                          - log(theta + mu_i) - (y_i + theta) / (theta + mu_i) ]
#   H(theta) = sum_i a_i [ psi'(y_i + theta) - psi'(theta) + 1 / theta
#                          - 2 / (theta + mu_i)
#                          + (y_i + theta) / (theta + mu_i)^2 ]
# psi = digamma(), psi' = trigamma().
#
# The code uses algebraically identical forms that avoid subtracting large,
# nearly equal numbers when theta is large:
#   log theta - log(theta + mu)                 = -log1p(mu / theta)
#   1 - (y + theta) / (theta + mu)              = (mu - y) / (theta + mu)
#   1/theta - 1/(theta + mu)                    = mu / (theta (theta + mu))
#   -1/(theta + mu) + (y + theta)/(theta + mu)^2 = (y - mu) / (theta + mu)^2
# so
#   U_i = psi(y + theta) - psi(theta) - log1p(mu / theta)
#         + (mu - y) / (theta + mu)
#   H_i = psi'(y + theta) - psi'(theta) + mu / (theta (theta + mu))
#         + (y - mu) / (theta + mu)^2
theta_score <- function(y, mu, theta, a) {
  U <- sum(a * (digamma(y + theta) - digamma(theta) - log1p(mu / theta) +
                  (mu - y) / (theta + mu)))
  H <- sum(a * (trigamma(y + theta) - trigamma(theta) +
                  mu / (theta * (theta + mu)) + (y - mu) / (theta + mu)^2))
  list(U = U, H = H)
}

# Newton-Raphson for theta with mu (i.e. beta) held fixed, Section 3.1(b).
# We update phi = log(theta), which keeps theta > 0. By the chain rule
# (d theta / d phi = theta):
#   U_phi = theta U,   H_phi = theta^2 H + theta U,
#   phi <- phi - U_phi / H_phi.
# Safeguards (not in the equations, needed for robustness):
#   * If H_phi >= 0 the Newton step would not point uphill, so we take a
#     unit step in the direction of the score instead.
#   * Step halving: shrink the step until the log-likelihood does not
#     decrease, so every update is an improvement.
#   * phi is capped at phi_max. As theta -> Inf the NB tends to the Poisson
#     and the log-likelihood flattens out, so theta could drift without
#     bound; nbreg() warns when this happens.
theta_newton <- function(y, mu, a, theta, tol = 1e-8, maxit = 100,
                         phi_max = log(1e10)) {
  phi <- log(theta)
  ll <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    s <- theta_score(y, mu, theta, a)
    U_phi <- theta * s$U
    H_phi <- theta^2 * s$H + theta * s$U
    step <- if (H_phi < 0) -U_phi / H_phi else sign(U_phi)
    # Step halving: accept the step once the log-likelihood is not lower.
    repeat {
      phi_new <- min(phi + step, phi_max)
      ll_new <- nb_loglik(y, mu, exp(phi_new), a)
      if (ll_new >= ll || abs(step) < tol) break
      step <- step / 2
    }
    phi <- phi_new
    theta <- exp(phi)
    ll <- ll_new
    if (abs(step) < tol) {
      converged <- TRUE
      break
    }
    if (phi >= phi_max) break  # hit the cap: theta is effectively infinite
  }
  list(theta = theta, iter = iter, converged = converged)
}

# IRLS for the zero part of the ZINB M-step, Section 3.2.
# Maximizes sum_i [ tau_i log pi_i + (1 - tau_i) log(1 - pi_i) ] over gamma,
# with logit pi_i = z_i' gamma. This is logistic regression with a
# fractional response tau_i in [0, 1]. For the logit link
# d pi / d eta = pi (1 - pi) and V(pi) = pi (1 - pi), so Section 3.2 gives
#   w_i = pi_i (1 - pi_i),
#   z_i = z_i' gamma + (tau_i - pi_i) / (pi_i (1 - pi_i)),
#   gamma(t+1) = (Z' W Z)^(-1) Z' W z.
# gamma = NULL starts from pi_i = (tau_i + 0.5) / 2, which keeps the start
# strictly inside (0, 1).
irls_logit <- function(Z, tau, gamma = NULL, tol = 1e-8, maxit = 100) {
  if (is.null(gamma)) {
    eta <- stats::qlogis((tau + 0.5) / 2)
  } else {
    eta <- drop(Z %*% gamma)
  }
  converged <- FALSE
  for (iter in seq_len(maxit)) {
    pi <- stats::plogis(eta)
    # w_i = pi_i (1 - pi_i), bounded away from 0 so the division below is
    # safe if some pi_i are numerically 0 or 1.
    w <- pmax(pi * (1 - pi), .Machine$double.eps)
    z <- eta + (tau - pi) / w
    gamma_new <- wls_solve(Z, z, w)
    eta <- drop(Z %*% gamma_new)
    if (!is.null(gamma) && max(abs(gamma_new - gamma)) < tol) {
      gamma <- gamma_new
      converged <- TRUE
      break
    }
    gamma <- gamma_new
  }
  list(gamma = gamma, pi = stats::plogis(eta), iter = iter,
       converged = converged)
}
