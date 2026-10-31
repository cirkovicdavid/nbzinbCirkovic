#' Zero-inflated negative binomial regression
#'
#' Fits a zero-inflated negative binomial (ZINB) model. With probability
#' \eqn{\pi_i}{pi_i} observation \eqn{i} is a structural zero; otherwise it
#' comes from an NB2 distribution with mean \eqn{\mu_i}{mu_i} and dispersion
#' \eqn{\theta}{theta}. The count model is
#' \eqn{\log\mu_i = x_i^\top\beta + o_i}{log mu_i = x_i' beta + o_i} and the
#' zero-inflation model is
#' \eqn{\mathrm{logit}\,\pi_i = z_i^\top\gamma}{logit pi_i = z_i' gamma}.
#'
#' The model is fitted by the EM algorithm. The E-step computes the
#' posterior probability \eqn{\tau_i}{tau_i} that each zero is a structural
#' zero. The M-step fits (i) a logistic regression of \eqn{\tau_i}{tau_i} on
#' \eqn{z_i} by IRLS and (ii) a weighted NB regression with prior weights
#' \eqn{1 - \tau_i}{1 - tau_i} by the same IRLS / Newton-Raphson algorithm as
#' [nbreg()]. Iterations stop when the observed-data log-likelihood changes by
#' less than `control$epsilon`. No general-purpose numerical maximizer is
#' used.
#'
#' EM does not produce standard errors. They are obtained by inverting the
#' negative Hessian of the observed-data log-likelihood with respect to
#' \eqn{(\beta, \gamma, \log\theta)}{(beta, gamma, log theta)} at the
#' estimates, computed with [numDeriv::hessian()].
#'
#' @param formula A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model
#'   covariates before `|`, zero-model covariates after it. Without `|`, the
#'   same covariates are used in both parts. An `offset()` term is allowed in
#'   the count part only.
#' @param data A data frame containing the variables in `formula`.
#' @param offset An optional offset for the count model (e.g. log exposure),
#'   evaluated in `data`.
#' @param control A list of settings from [zinbreg_control()].
#'
#' @return An object of class `"zinbreg"`, a list with components
#' \describe{
#'   \item{coefficients}{a list with elements `count` (\eqn{\hat\beta}{beta-hat})
#'     and `zero` (\eqn{\hat\gamma}{gamma-hat})}
#'   \item{vcov}{variance-covariance matrix of
#'     \eqn{(\hat\beta, \hat\gamma)}{(beta-hat, gamma-hat)}; names are prefixed
#'     `count_` and `zero_`}
#'   \item{theta, SE.theta, SE.logtheta}{estimated dispersion parameter and
#'     the standard errors of \eqn{\hat\theta}{theta-hat} and
#'     \eqn{\log\hat\theta}{log theta-hat}}
#'   \item{vcov_full}{variance-covariance matrix of
#'     \eqn{(\hat\beta, \hat\gamma, \log\hat\theta)}{(beta-hat, gamma-hat, log theta-hat)}}
#'   \item{loglik, loglik_path}{maximized log-likelihood and its value at
#'     every EM iteration}
#'   \item{iter, converged}{number of EM iterations and convergence flag}
#'   \item{fitted.values}{\eqn{(1 - \hat\pi_i)\hat\mu_i}{(1 - pi-hat) mu-hat},
#'     the fitted mean of \eqn{Y_i}}
#'   \item{mu, pi, tau}{fitted count means, structural-zero probabilities, and
#'     final E-step posterior probabilities}
#'   \item{y, offset, nobs, df.residual}{response, offset, number of
#'     observations, and residual degrees of freedom}
#'   \item{call, formula, terms, model, xlevels, contrasts, control}{information
#'     used to print the fit and to predict for new data}
#' }
#'
#' @seealso [zinbreg_control()], [nbreg()]
#'
#' @examples
#' # Simulated ZINB data
#' set.seed(1)
#' n <- 500
#' x <- rnorm(n)
#' z <- rnorm(n)
#' structural_zero <- rbinom(n, 1, plogis(-1 + 0.8 * z))
#' y <- ifelse(structural_zero == 1, 0,
#'             rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)))
#' fit <- zinbreg(y ~ x | z, data = data.frame(x, z, y))
#' fit$coefficients
#' fit$theta
#'
#' # Publication counts of biochemistry PhD students
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   fit2 <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
#'                   data = bioChemists)
#'   fit2$coefficients
#' }
#'
#' @importFrom Formula as.Formula
#' @importFrom stats terms
#' @export
zinbreg <- function(formula, data, offset, control = zinbreg_control()) {
  cl <- match.call()

  # Two-part formula. Without "|", reuse the count covariates for the zero
  # part: y ~ x becomes y ~ x | x.
  ff <- as.Formula(formula)
  if (length(ff)[2L] == 1L) {
    ff <- as.Formula(formula(ff), formula(ff, lhs = 0L, rhs = 1L))
  } else if (length(ff)[2L] > 2L) {
    stop("the formula must have at most two parts: y ~ x | z", call. = FALSE)
  }

  # Model frame (same idiom as nbreg(), but with the two-part formula;
  # model.frame() dispatches to the Formula package's method).
  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf$formula <- ff
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Separate terms and model matrices for the count part (rhs = 1) and the
  # zero part (rhs = 2).
  mtX <- terms(ff, data = mf, rhs = 1L)
  mtZ <- terms(ff, data = mf, rhs = 2L)
  if (!is.null(attr(mtZ, "offset"))) {
    stop("offsets are only supported in the count part of the formula",
         call. = FALSE)
  }
  X <- model.matrix(mtX, mf)
  Z <- model.matrix(mtZ, mf)
  y <- model.response(mf)
  off <- model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(X))

  check_counts(y)
  if (!any(y == 0)) {
    stop("the response has no zeros, so a zero-inflated model is not ",
         "identifiable", call. = FALSE)
  }
  if (ncol(X) == 0L || ncol(Z) == 0L) {
    stop("both the count and the zero part need at least one coefficient",
         call. = FALSE)
  }
  if (nrow(X) <= ncol(X) + ncol(Z) + 1L) {
    stop("more parameters than observations", call. = FALSE)
  }

  fit <- zinb_fit(X, Z, y, off, control)

  if (!fit$converged) {
    warning("zinbreg() did not converge in ", control$maxit,
            " EM iterations; consider increasing zinbreg_control(maxit = )",
            call. = FALSE)
  }
  if (fit$decreased) {
    warning("the log-likelihood decreased during EM; ",
            "check the data and starting values", call. = FALSE)
  }
  if (fit$theta > control$theta_warn) {
    warning("theta estimate is very large (", signif(fit$theta, 3), "): ",
            "the count part shows little overdispersion, so a ",
            "zero-inflated Poisson model may be adequate", call. = FALSE)
  }

  # Standard errors (Section 3.2): invert the observed information, i.e.
  # minus the Hessian of the observed-data log-likelihood with respect to
  # par = (beta, gamma, log theta), at the estimates.
  p <- ncol(X)
  q <- ncol(Z)
  par_hat <- c(fit$beta, fit$gamma, log(fit$theta))
  par_names <- c(paste0("count_", colnames(X)), paste0("zero_", colnames(Z)),
                 "Log(theta)")
  H <- numDeriv::hessian(zinb_loglik_par, par_hat, X = X, Z = Z, y = y,
                         offset = off)
  vc_full <- tryCatch(solve(-H), error = function(e) {
    warning("the Hessian is singular; standard errors are not available",
            call. = FALSE)
    matrix(NA_real_, p + q + 1L, p + q + 1L)
  })
  dimnames(vc_full) <- list(par_names, par_names)
  idx <- seq_len(p + q)
  se_logtheta <- sqrt(vc_full[p + q + 1L, p + q + 1L])

  structure(list(
    coefficients = list(count = fit$beta, zero = fit$gamma),
    vcov = vc_full[idx, idx, drop = FALSE],
    vcov_full = vc_full,
    theta = fit$theta,
    SE.logtheta = se_logtheta,
    # Delta method: SE(theta) = theta * SE(log theta)
    SE.theta = fit$theta * se_logtheta,
    loglik = fit$loglik,
    loglik_path = fit$loglik_path,
    iter = fit$iter,
    converged = fit$converged,
    fitted.values = stats::setNames((1 - fit$pi) * fit$mu, rownames(mf)),
    mu = stats::setNames(fit$mu, rownames(mf)),
    pi = stats::setNames(fit$pi, rownames(mf)),
    tau = stats::setNames(fit$tau, rownames(mf)),
    y = y,
    offset = off,
    nobs = length(y),
    df.residual = length(y) - (p + q + 1L),
    call = cl,
    formula = ff,
    terms = list(count = mtX, zero = mtZ, full = attr(mf, "terms")),
    model = mf,
    xlevels = list(count = .getXlevels(mtX, mf), zero = .getXlevels(mtZ, mf)),
    contrasts = list(count = attr(X, "contrasts"),
                     zero = attr(Z, "contrasts")),
    control = control
  ), class = "zinbreg")
}

#' Control settings for zinbreg()
#'
#' @param epsilon Convergence tolerance. EM stops when the absolute change in
#'   the observed-data log-likelihood is below `epsilon`. The same tolerance
#'   is used inside the M-steps.
#' @param maxit Maximum number of EM iterations.
#' @param irls_maxit Maximum number of IRLS iterations in each M-step fit.
#' @param theta_maxit Maximum number of Newton-Raphson iterations for theta
#'   in each count M-step.
#' @param nb_maxit Maximum number of outer (theta / beta) iterations in each
#'   count M-step.
#' @param theta_warn Warn if the estimated theta exceeds this value.
#' @param trace Logical; if `TRUE`, print the observed-data log-likelihood
#'   and its change at each EM iteration.
#'
#' @return A list of control settings.
#' @examples
#' zinbreg_control(epsilon = 1e-10, trace = TRUE)
#' @export
zinbreg_control <- function(epsilon = 1e-8, maxit = 5000, irls_maxit = 100,
                            theta_maxit = 100, nb_maxit = 100,
                            theta_warn = 1e4, trace = FALSE) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("'epsilon' must be a single positive number", call. = FALSE)
  }
  for (nm in c("maxit", "irls_maxit", "theta_maxit", "nb_maxit")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || v < 1) {
      stop("'", nm, "' must be a single number >= 1", call. = FALSE)
    }
  }
  list(epsilon = epsilon, maxit = as.integer(maxit),
       irls_maxit = as.integer(irls_maxit),
       theta_maxit = as.integer(theta_maxit),
       nb_maxit = as.integer(nb_maxit), theta_warn = theta_warn,
       trace = isTRUE(trace))
}

# EM algorithm for the ZINB model, Section 3.2.
zinb_fit <- function(X, Z, y, offset, control) {
  tol <- control$epsilon
  is0 <- y == 0
  # Settings for the count M-step (an NB fit with prior weights 1 - tau).
  nb_control <- list(epsilon = tol, maxit = control$nb_maxit,
                     irls_maxit = control$irls_maxit,
                     theta_maxit = control$theta_maxit, trace = FALSE)

  # Starting values, Section 3.2:
  #   (beta, theta): an ordinary NB fit (all prior weights 1);
  #   gamma: logistic IRLS of I(y_i = 0) on z_i.
  nb0 <- nb_fit(X, y, offset, rep(1, length(y)), nb_control)
  beta <- nb0$beta
  theta <- nb0$theta
  mu <- nb0$mu
  lg0 <- irls_logit(Z, as.numeric(is0), tol = tol,
                    maxit = control$irls_maxit)
  gamma <- lg0$gamma
  pi <- lg0$pi

  ll <- zinb_loglik(y, mu, theta, pi)
  ll_path <- ll
  converged <- FALSE
  decreased <- FALSE
  for (iter in seq_len(control$maxit)) {
    # E-step: posterior probability of a structural zero,
    #   tau_i = pi_i / (pi_i + (1 - pi_i) f_NB(0; mu_i, theta))  if y_i = 0,
    #   tau_i = 0                                                 if y_i > 0
    # (a positive count cannot be a structural zero).
    f0 <- stats::dnbinom(0, size = theta, mu = mu)
    tau <- ifelse(is0, pi / (pi + (1 - pi) * f0), 0)

    # M-step, zero part: maximize sum tau log pi + (1 - tau) log(1 - pi)
    # over gamma (logistic IRLS with fractional response tau).
    lg <- irls_logit(Z, tau, gamma = gamma, tol = tol,
                     maxit = control$irls_maxit)
    gamma <- lg$gamma
    pi <- lg$pi

    # M-step, count part: maximize sum (1 - tau) log f_NB(y; mu, theta) over
    # (beta, theta) = the NB fit of Section 3.1 with prior weights
    # a_i = 1 - tau_i, warm-started at the current (beta, theta).
    nb <- nb_fit(X, y, offset, a = 1 - tau, control = nb_control,
                 start = list(beta = beta, theta = theta))
    beta <- nb$beta
    theta <- nb$theta
    mu <- nb$mu

    # Convergence: change in the observed-data log-likelihood. EM must
    # never decrease it; a decrease beyond rounding signals a bug.
    ll_new <- zinb_loglik(y, mu, theta, pi)
    if (ll_new < ll - 1e-6) decreased <- TRUE
    if (control$trace) {
      message(sprintf("zinbreg EM iter %4d: logLik = %.10f  change = %.3e",
                      iter, ll_new, ll_new - ll))
    }
    ll_path <- c(ll_path, ll_new)
    done <- abs(ll_new - ll) < tol
    ll <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }
  # Final E-step so tau matches the reported estimates.
  f0 <- stats::dnbinom(0, size = theta, mu = mu)
  tau <- ifelse(is0, pi / (pi + (1 - pi) * f0), 0)

  list(beta = beta, gamma = gamma, theta = theta, mu = mu, pi = pi,
       tau = tau, loglik = ll, loglik_path = ll_path, iter = iter,
       converged = converged, decreased = decreased)
}
