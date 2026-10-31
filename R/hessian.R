# Analytic score and Hessian of the ZINB observed-data log-likelihood
# (Section 2.2) with respect to par = (beta, gamma, phi), phi = log(theta).
# Used for the standard errors in Section 3.2 ("analytic Hessian").
#
# Notation, per observation i:
#   eta_i  = x_i' beta + o_i,  mu_i = exp(eta_i)        (count part)
#   zeta_i = z_i' gamma,       pi_i = plogis(zeta_i)    (zero part)
#   phi    = log(theta)
#
# Chain rule. l_i depends on beta only through eta_i and on gamma only
# through zeta_i, with d eta_i / d beta = x_i and d zeta_i / d gamma = z_i.
# So, writing l_ab for the scalar second derivative of l_i in (a, b):
#   d^2 l / d beta d beta'   = sum_i l_{eta eta,i}   x_i x_i' = X' diag(.) X
#   d^2 l / d beta d gamma'  = sum_i l_{eta zeta,i}  x_i z_i' = X' diag(.) Z
#   d^2 l / d beta d phi     = sum_i l_{eta phi,i}   x_i      = X' (.)
#   d^2 l / d gamma d gamma' = sum_i l_{zeta zeta,i} z_i z_i' = Z' diag(.) Z
#   d^2 l / d gamma d phi    = sum_i l_{zeta phi,i}  z_i      = Z' (.)
#   d^2 l / d phi^2          = sum_i l_{phi phi,i}
# The functions below compute the six scalar derivatives for every i.
#
# ---------------------------------------------------------------------------
# Case y_i > 0:  l_i = log(1 - pi_i) + log f_NB(y_i; mu_i, theta)
#   The two terms depend on different parameters, so every cross derivative
#   between zeta and (eta, phi) is 0.
#   Zero part: log(1 - pi) = -log(1 + e^zeta)
#     l_zeta = -pi,   l_zeta zeta = -pi (1 - pi)
#   Count part: log f_NB = ... + y log mu - (theta + y) log(theta + mu), so
#     l_eta     = y - (theta + y) mu / (theta + mu) = theta (y - mu) / (theta + mu)
#     l_eta eta = -theta mu (theta + y) / (theta + mu)^2
#     l_eta theta = d/d theta [theta (y - mu) / (theta + mu)]
#                 = mu (y - mu) / (theta + mu)^2
#     l_eta phi = theta * l_eta theta   (d theta / d phi = theta)
#     l_phi     = theta U_i,   l_phi phi = theta^2 H_i + theta U_i,
#   with U_i, H_i the per-observation theta score and second derivative of
#   Section 3.1(b) (same formulas as theta_score(), without the sum).
#
# ---------------------------------------------------------------------------
# Case y_i = 0:  l_i = log L_i,  L_i = pi_i + (1 - pi_i) f0_i,
#   f0 = f_NB(0; mu, theta) = (theta / (theta + mu))^theta,
#   log f0 = -theta log(1 + mu / theta).
# For any parameters a, b:  l_a = L_a / L,  l_ab = L_ab / L - L_a L_b / L^2.
#
# Derivatives of g = log f0 (then f0_a = f0 g_a, f0_ab = f0 (g_a g_b + g_ab)):
#   g_eta     = -theta mu / (theta + mu)
#   g_theta   = -log(1 + mu / theta) + mu / (theta + mu)
#   g_eta eta = -theta^2 mu / (theta + mu)^2
#   g_eta theta = -mu^2 / (theta + mu)^2
#   g_theta theta = mu^2 / (theta (theta + mu)^2)
#   On the phi scale: g_phi = theta g_theta,  g_eta phi = theta g_eta theta,
#   g_phi phi = theta^2 g_theta theta + theta g_theta.
#
# Derivatives of L. With s = d pi / d zeta = pi (1 - pi) and
# d^2 pi / d zeta^2 = s (1 - 2 pi):
#   L_zeta      = s (1 - f0)
#   L_zeta zeta = s (1 - 2 pi) (1 - f0)
#   L_a         = (1 - pi) f0 g_a                  (a = eta or phi)
#   L_ab        = (1 - pi) f0 (g_a g_b + g_ab)     (a, b in {eta, phi})
#   L_zeta a    = -s f0 g_a
# ---------------------------------------------------------------------------

# Per-observation first and second derivatives of l_i with respect to
# (eta_i, zeta_i, phi). Returns a list of n-vectors.
zinb_obs_derivs <- function(y, mu, theta, pi) {
  is0 <- y == 0
  tm <- theta + mu

  # ---- y > 0 ---------------------------------------------------------------
  # Per-observation theta score U_i and second derivative H_i (Section 3.1(b),
  # in the cancellation-free form explained in theta_score()).
  U_i <- digamma(y + theta) - digamma(theta) - log1p(mu / theta) +
    (mu - y) / tm
  H_i <- trigamma(y + theta) - trigamma(theta) + mu / (theta * tm) +
    (y - mu) / tm^2
  pos <- list(
    eta = theta * (y - mu) / tm,
    zeta = -pi,
    phi = theta * U_i,
    eta_eta = -theta * mu * (theta + y) / tm^2,
    eta_zeta = rep(0, length(y)),
    eta_phi = theta * mu * (y - mu) / tm^2,
    zeta_zeta = -pi * (1 - pi),
    zeta_phi = rep(0, length(y)),
    phi_phi = theta^2 * H_i + theta * U_i
  )

  # ---- y = 0 ---------------------------------------------------------------
  f0 <- stats::dnbinom(0, size = theta, mu = mu)
  s <- pi * (1 - pi)
  # derivatives of g = log f0
  g_eta <- -theta * mu / tm
  g_theta <- -log1p(mu / theta) + mu / tm
  g_phi <- theta * g_theta
  g_eta_eta <- -theta^2 * mu / tm^2
  g_eta_phi <- theta * (-mu^2 / tm^2)
  g_phi_phi <- theta^2 * (mu^2 / (theta * tm^2)) + theta * g_theta
  # derivatives of L = pi + (1 - pi) f0
  L <- pi + (1 - pi) * f0
  c1 <- (1 - pi) * f0
  L_eta <- c1 * g_eta
  L_phi <- c1 * g_phi
  L_zeta <- s * (1 - f0)
  L_eta_eta <- c1 * (g_eta^2 + g_eta_eta)
  L_eta_phi <- c1 * (g_eta * g_phi + g_eta_phi)
  L_phi_phi <- c1 * (g_phi^2 + g_phi_phi)
  L_zeta_zeta <- s * (1 - 2 * pi) * (1 - f0)
  L_eta_zeta <- -s * f0 * g_eta
  L_zeta_phi <- -s * f0 * g_phi
  # l = log L:  l_a = L_a / L,  l_ab = L_ab / L - L_a L_b / L^2
  zero <- list(
    eta = L_eta / L,
    zeta = L_zeta / L,
    phi = L_phi / L,
    eta_eta = L_eta_eta / L - L_eta^2 / L^2,
    eta_zeta = L_eta_zeta / L - L_eta * L_zeta / L^2,
    eta_phi = L_eta_phi / L - L_eta * L_phi / L^2,
    zeta_zeta = L_zeta_zeta / L - L_zeta^2 / L^2,
    zeta_phi = L_zeta_phi / L - L_zeta * L_phi / L^2,
    phi_phi = L_phi_phi / L - L_phi^2 / L^2
  )

  # Pick the y = 0 or y > 0 formula for each observation.
  Map(function(a, b) ifelse(is0, a, b), zero, pos)
}

# Unpack par = (beta, gamma, log theta) into mu, pi, theta.
zinb_unpack <- function(par, X, Z, offset) {
  p <- ncol(X)
  q <- ncol(Z)
  theta <- exp(par[p + q + 1])
  list(mu = exp(drop(X %*% par[seq_len(p)]) + offset),
       pi = stats::plogis(drop(Z %*% par[p + seq_len(q)])),
       theta = theta)
}

# Score vector d l / d par, par = (beta, gamma, log theta).
zinb_score <- function(par, X, Z, y, offset) {
  u <- zinb_unpack(par, X, Z, offset)
  d <- zinb_obs_derivs(y, u$mu, u$theta, u$pi)
  c(crossprod(X, d$eta), crossprod(Z, d$zeta), sum(d$phi))
}

# Hessian d^2 l / d par d par', assembled block by block from the chain-rule
# formulas at the top of this file.
zinb_hessian <- function(par, X, Z, y, offset) {
  u <- zinb_unpack(par, X, Z, offset)
  d <- zinb_obs_derivs(y, u$mu, u$theta, u$pi)
  H_bb <- crossprod(X, d$eta_eta * X)    # X' diag(l_eta eta) X
  H_bg <- crossprod(X, d$eta_zeta * Z)   # X' diag(l_eta zeta) Z
  H_bp <- crossprod(X, d$eta_phi)        # X' l_eta phi
  H_gg <- crossprod(Z, d$zeta_zeta * Z)  # Z' diag(l_zeta zeta) Z
  H_gp <- crossprod(Z, d$zeta_phi)       # Z' l_zeta phi
  H_pp <- sum(d$phi_phi)
  rbind(cbind(H_bb, H_bg, H_bp),
        cbind(t(H_bg), H_gg, H_gp),
        cbind(t(H_bp), t(H_gp), H_pp))
}
