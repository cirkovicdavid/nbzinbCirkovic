#' Negative binomial regression
#'
#' Fits a negative binomial (NB2) regression model with log link,
#' \eqn{\log \mu_i = x_i^\top \beta + o_i}{log mu_i = x_i' beta + o_i} and
#' \eqn{\mathrm{Var}(Y_i) = \mu_i + \mu_i^2/\theta}{Var(Y_i) = mu_i + mu_i^2 / theta}.
#' The coefficients \eqn{\beta}{beta} are estimated by iteratively reweighted
#' least squares (IRLS) with \eqn{\theta}{theta} held fixed, alternated with
#' Newton-Raphson updates of \eqn{\log\theta}{log(theta)} with
#' \eqn{\beta}{beta} held fixed, until the log-likelihood stops changing.
#' No general-purpose numerical maximizer is used.
#'
#' The algorithm starts from a Poisson fit (IRLS with
#' \eqn{\theta = \infty}{theta = Inf}) and the moment estimate
#' \eqn{\theta_0 = n / \sum_i (y_i/\mu_i - 1)^2}{theta_0 = n / sum((y / mu - 1)^2)}.
#' Standard errors of \eqn{\hat\beta}{beta-hat} come from
#' \eqn{(X^\top \hat W X)^{-1}}{(X' W X)^(-1)}, treating \eqn{\theta}{theta}
#' as known, and the standard error of \eqn{\hat\theta}{theta-hat} is
#' \eqn{\{-H(\hat\theta)\}^{-1/2}}{(-H(theta-hat))^(-1/2)}, where \eqn{H} is
#' the second derivative of the log-likelihood in \eqn{\theta}{theta}.
#'
#' @param formula A model formula, `response ~ covariates`. Offsets may be
#'   given in the formula with `offset()`.
#' @param data A data frame containing the variables in `formula`.
#' @param offset An optional offset vector (e.g. log exposure), evaluated in
#'   `data`. Added to any `offset()` terms in the formula.
#' @param control A list of settings from [nbreg_control()].
#'
#' @return An object of class `"nbreg"`, a list with components
#' \describe{
#'   \item{coefficients}{estimated regression coefficients \eqn{\hat\beta}{beta-hat}}
#'   \item{vcov}{variance-covariance matrix of \eqn{\hat\beta}{beta-hat}}
#'   \item{theta, SE.theta}{estimated dispersion parameter and its standard error}
#'   \item{loglik}{maximized log-likelihood}
#'   \item{iter, converged}{number of outer iterations and whether the
#'     log-likelihood converged}
#'   \item{fitted.values, linear.predictors}{\eqn{\hat\mu_i}{mu-hat} and
#'     \eqn{\hat\eta_i}{eta-hat}}
#'   \item{y, offset, nobs, df.residual}{response, offset, number of
#'     observations, and residual degrees of freedom}
#'   \item{call, terms, model, xlevels, contrasts, control}{information used
#'     to print the fit and to predict for new data}
#' }
#'
#' @seealso [nbreg_control()], [zinbreg()]
#'
#' @examples
#' # Simulated data with known beta = (0.5, 0.3) and theta = 2
#' set.seed(1)
#' n <- 300
#' x <- rnorm(n)
#' y <- rnbinom(n, size = 2, mu = exp(0.5 + 0.3 * x))
#' fit <- nbreg(y ~ x, data = data.frame(x, y))
#' fit$coefficients
#' fit$theta
#'
#' # Publication counts of biochemistry PhD students
#' if (requireNamespace("pscl", quietly = TRUE)) {
#'   data("bioChemists", package = "pscl")
#'   fit2 <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
#'   cbind(estimate = fit2$coefficients, se = sqrt(diag(fit2$vcov)))
#' }
#'
#' @importFrom stats model.frame model.matrix model.response model.offset
#'   .getXlevels
#' @export
nbreg <- function(formula, data, offset, control = nbreg_control()) {
  cl <- match.call()

  # Build the model frame from formula, data and offset (the standard idiom
  # used by R's modelling functions: keep only those arguments of the call
  # and evaluate them with model.frame()).
  mf <- match.call(expand.dots = FALSE)
  keep <- match(c("formula", "data", "offset"), names(mf), 0L)
  mf <- mf[c(1L, keep)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  y <- model.response(mf)
  X <- model.matrix(mt, mf)
  # model.offset() adds the offset argument and any offset() terms.
  off <- model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(X))

  check_counts(y)
  if (ncol(X) == 0L) stop("the model must have at least one coefficient",
                          call. = FALSE)
  if (nrow(X) <= ncol(X)) stop("more coefficients than observations",
                               call. = FALSE)

  a <- rep(1, length(y))  # ordinary fit: all prior weights a_i = 1
  fit <- nb_fit(X, y, off, a, control)

  if (!fit$converged) {
    warning("nbreg() did not converge in ", control$maxit, " iterations; ",
            "consider increasing nbreg_control(maxit = )", call. = FALSE)
  }
  if (fit$theta > control$theta_warn) {
    warning("theta estimate is very large (", signif(fit$theta, 3), "): ",
            "the data show little overdispersion, so a Poisson model ",
            "may be adequate", call. = FALSE)
  }

  # Standard errors, Section 3.1(d).
  # Var(beta-hat) = (X' W X)^(-1) with W evaluated at (beta-hat, theta-hat),
  # treating theta as known.
  w <- irls_weights(fit$mu, fit$theta, a)
  vc <- chol2inv(chol(crossprod(X, w * X)))
  dimnames(vc) <- list(colnames(X), colnames(X))
  # SE(theta-hat) = (-H(theta-hat))^(-1/2). At a proper maximum H < 0. In the
  # Poisson limit (theta -> Inf) H is ~0 and rounding can make it >= 0; theta
  # is then not estimable and its SE is reported as NA.
  H <- theta_score(y, fit$mu, fit$theta, a)$H
  se_theta <- if (H < 0) 1 / sqrt(-H) else NA_real_

  structure(list(
    coefficients = fit$beta,
    vcov = vc,
    theta = fit$theta,
    SE.theta = se_theta,
    loglik = fit$loglik,
    iter = fit$iter,
    converged = fit$converged,
    fitted.values = stats::setNames(fit$mu, rownames(mf)),
    linear.predictors = stats::setNames(fit$eta, rownames(mf)),
    y = y,
    offset = off,
    nobs = length(y),
    df.residual = length(y) - ncol(X),
    call = cl,
    terms = mt,
    model = mf,
    xlevels = .getXlevels(mt, mf),
    contrasts = attr(X, "contrasts"),
    control = control
  ), class = "nbreg")
}

#' Control settings for nbreg()
#'
#' @param epsilon Convergence tolerance. The outer loop stops when the
#'   absolute change in the log-likelihood is below `epsilon`; the inner
#'   IRLS and Newton-Raphson loops stop when their parameter updates are
#'   below `epsilon`.
#' @param maxit Maximum number of outer iterations (each one is a
#'   Newton-Raphson update of theta followed by an IRLS fit of beta).
#' @param irls_maxit Maximum number of IRLS iterations per outer iteration.
#' @param theta_maxit Maximum number of Newton-Raphson iterations per outer
#'   iteration.
#' @param theta_warn Warn if the estimated theta exceeds this value
#'   (the Poisson limit).
#' @param trace Logical; if `TRUE`, print the log-likelihood and theta at
#'   each outer iteration.
#'
#' @return A list of control settings.
#' @examples
#' nbreg_control(epsilon = 1e-10, trace = TRUE)
#' @export
nbreg_control <- function(epsilon = 1e-8, maxit = 100, irls_maxit = 100,
                          theta_maxit = 100, theta_warn = 1e4,
                          trace = FALSE) {
  if (!is.numeric(epsilon) || length(epsilon) != 1L || epsilon <= 0) {
    stop("'epsilon' must be a single positive number", call. = FALSE)
  }
  for (nm in c("maxit", "irls_maxit", "theta_maxit")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || v < 1) {
      stop("'", nm, "' must be a single number >= 1", call. = FALSE)
    }
  }
  list(epsilon = epsilon, maxit = as.integer(maxit),
       irls_maxit = as.integer(irls_maxit),
       theta_maxit = as.integer(theta_maxit), theta_warn = theta_warn,
       trace = isTRUE(trace))
}

# NB fitting engine, Section 3.1(c). Shared by nbreg() and by the count
# M-step of zinbreg() (with a_i = 1 - tau_i). Returns estimates only; the
# caller computes standard errors and issues warnings.
#   start  optional list(beta, theta) for a warm start (used inside EM)
nb_fit <- function(X, y, offset, a, control, start = NULL) {
  tol <- control$epsilon
  if (is.null(start)) {
    # Start from a Poisson fit: IRLS with theta = Inf (w_i = a_i mu_i).
    pois <- irls_nb(X, y, offset, a, theta = Inf, beta = NULL, tol = tol,
                    maxit = control$irls_maxit)
    beta <- pois$beta
    mu <- pois$mu
    # Moment estimate of theta: E[(y - mu)^2 / mu^2] = 1/mu + 1/theta,
    # approximated as theta_0 = sum(a) / sum(a (y / mu - 1)^2).
    theta <- sum(a) / sum(a * (y / mu - 1)^2)
  } else {
    beta <- start$beta
    theta <- start$theta
    mu <- exp(drop(X %*% beta) + offset)
  }
  eta <- log(mu)

  ll <- nb_loglik(y, mu, theta, a)
  converged <- FALSE
  for (iter in seq_len(control$maxit)) {
    # (b) Newton-Raphson for log(theta), beta fixed
    theta <- theta_newton(y, mu, a, theta, tol = tol,
                          maxit = control$theta_maxit)$theta
    # (a) IRLS for beta, theta fixed (warm-started at the current beta)
    ir <- irls_nb(X, y, offset, a, theta, beta = beta, tol = tol,
                  maxit = control$irls_maxit)
    beta <- ir$beta
    mu <- ir$mu
    eta <- ir$eta
    # (c) Stop when the log-likelihood changes by less than epsilon
    ll_new <- nb_loglik(y, mu, theta, a)
    if (control$trace) {
      message(sprintf("nbreg iter %3d: logLik = %.10f  theta = %.8g",
                      iter, ll_new, theta))
    }
    done <- abs(ll_new - ll) < tol
    ll <- ll_new
    if (done) {
      converged <- TRUE
      break
    }
  }
  list(beta = beta, theta = theta, mu = mu, eta = eta, loglik = ll,
       iter = iter, converged = converged)
}

# Input check: the response must be non-negative whole numbers.
check_counts <- function(y) {
  if (is.null(y)) stop("the formula needs a response", call. = FALSE)
  if (!is.numeric(y) || !is.null(dim(y))) {
    stop("the response must be a numeric vector of counts", call. = FALSE)
  }
  if (any(y < 0) || any(abs(y - round(y)) > 1e-8)) {
    stop("the response must contain non-negative integer counts",
         call. = FALSE)
  }
  invisible(TRUE)
}
