# S3 methods for "nbreg" objects (Section 4.2).

#' Methods for nbreg fits
#'
#' Standard methods for objects returned by [nbreg()].
#'
#' @param object,x An object of class `"nbreg"` (or `"summary.nbreg"` for
#'   the `print` method of a summary).
#' @param type For `residuals()`: `"pearson"` (default),
#'   \eqn{(y_i - \hat\mu_i) / \sqrt{\hat\mu_i + \hat\mu_i^2/\hat\theta}}{(y - mu) / sqrt(mu + mu^2 / theta)},
#'   or `"response"`, \eqn{y_i - \hat\mu_i}{y - mu}. For `predict()`:
#'   `"response"` (default) for the mean \eqn{\mu}{mu}, or `"link"` for the
#'   linear predictor \eqn{\eta = \log\mu}{eta = log(mu)}.
#' @param newdata Optional data frame of covariate values for prediction.
#'   If omitted, predictions for the data used in the fit are returned.
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (currently ignored, except by `print()`
#'   methods which pass them to [stats::printCoefmat()]).
#'
#' @return `coef()` the coefficient vector; `vcov()` its covariance matrix;
#'   `logLik()` a `"logLik"` object with `df` = number of coefficients + 1
#'   (for \eqn{\theta}{theta}); `fitted()`, `residuals()` and `predict()`
#'   numeric vectors; `nobs()` the number of observations; `summary()` an
#'   object of class `"summary.nbreg"` with the coefficient table.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200))
#' d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.3 * d$x))
#' fit <- nbreg(y ~ x, data = d)
#' fit
#' summary(fit)
#' AIC(fit)
#' head(residuals(fit, type = "pearson"))
#' predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#' @name nbreg-methods
NULL

#' @rdname nbreg-methods
#' @export
print.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Negative binomial regression (IRLS + Newton-Raphson for theta)\n\n")
  cat("Coefficients:\n")
  print(format(x$coefficients, digits = digits), quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Log-likelihood:", format(x$loglik, digits = digits + 2),
      "on", length(x$coefficients) + 1L, "df\n")
  if (!x$converged) cat("Warning: the algorithm did not converge\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
summary.nbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov))
  structure(list(
    call = object$call,
    coefficients = coef_table(object$coefficients, se),
    theta = object$theta,
    SE.theta = object$SE.theta,
    loglik = object$loglik,
    df = length(object$coefficients) + 1L,
    aic = stats::AIC(object),
    iter = object$iter,
    converged = object$converged,
    residuals = stats::residuals(object, type = "pearson")
  ), class = "summary.nbreg")
}

#' @rdname nbreg-methods
#' @export
print.summary.nbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                                ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  print(stats::setNames(stats::quantile(x$residuals),
                        c("Min", "1Q", "Median", "3Q", "Max")),
        digits = digits)
  cat("\nCoefficients (negative binomial with log link):\n")
  stats::printCoefmat(x$coefficients, digits = digits, ...)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Std. Error:", format(x$SE.theta, digits = digits), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2), "on",
      x$df, "df   AIC:", format(x$aic, digits = digits + 2), "\n")
  cat("Iterations:", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n")
  invisible(x)
}

#' @rdname nbreg-methods
#' @export
coef.nbreg <- function(object, ...) object$coefficients

#' @rdname nbreg-methods
#' @export
vcov.nbreg <- function(object, ...) object$vcov

#' @rdname nbreg-methods
#' @export
logLik.nbreg <- function(object, ...) {
  make_logLik(object$loglik, df = length(object$coefficients) + 1L,
              nobs = object$nobs)
}

#' @rdname nbreg-methods
#' @export
nobs.nbreg <- function(object, ...) object$nobs

#' @rdname nbreg-methods
#' @export
fitted.nbreg <- function(object, ...) object$fitted.values

#' @rdname nbreg-methods
#' @export
residuals.nbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  mu <- object$fitted.values
  r <- object$y - mu
  if (type == "pearson") {
    # Var(Y) = mu + mu^2 / theta (Section 2.1)
    r <- r / sqrt(mu + mu^2 / object$theta)
  }
  r
}

#' @rdname nbreg-methods
#' @export
predict.nbreg <- function(object, newdata = NULL,
                          type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    eta <- object$linear.predictors
  } else {
    des <- new_design(object$terms, newdata, object$xlevels,
                      object$contrasts, object$call)
    # eta = x' beta + o
    eta <- drop(des$X %*% object$coefficients) + des$offset
    names(eta) <- rownames(newdata)
  }
  if (type == "response") exp(eta) else eta
}
