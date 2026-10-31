# S3 methods for "zinbreg" objects (Section 4.2).

#' Methods for zinbreg fits
#'
#' Standard methods for objects returned by [zinbreg()].
#'
#' @param object,x An object of class `"zinbreg"` (or `"summary.zinbreg"`
#'   for the `print` method of a summary).
#' @param model For `coef()` and `vcov()`: `"full"` (default; count and
#'   zero coefficients, named `count_*` and `zero_*`), `"count"`, or
#'   `"zero"`.
#' @param type For `residuals()`: `"pearson"` (default) or `"response"`.
#'   The response residual is \eqn{y_i - (1 - \hat\pi_i)\hat\mu_i}{y - (1 - pi) mu};
#'   the Pearson residual divides it by
#'   \eqn{\sqrt{(1-\pi_i)\mu_i(1 + \mu_i(\pi_i + 1/\theta))}}{sqrt((1 - pi) mu (1 + mu (pi + 1/theta)))},
#'   the standard deviation of a ZINB variable. For `predict()`:
#'   `"response"` (default) for the mean \eqn{(1 - \pi)\mu}{(1 - pi) mu},
#'   `"count"` for the count-model mean \eqn{\mu}{mu}, or `"zero"` for the
#'   structural-zero probability \eqn{\pi}{pi}.
#' @param newdata Optional data frame of covariate values for prediction.
#' @param digits Number of significant digits to print.
#' @param ... Further arguments (currently ignored, except by `print()`
#'   methods which pass them to [stats::printCoefmat()]).
#'
#' @return `coef()` and `vcov()` the coefficients and their covariance
#'   matrix; `logLik()` a `"logLik"` object with `df` = number of count and
#'   zero coefficients + 1 (for \eqn{\theta}{theta}); `fitted()`,
#'   `residuals()` and `predict()` numeric vectors; `nobs()` the number of
#'   observations; `summary()` an object of class `"summary.zinbreg"`.
#'
#' @examples
#' set.seed(1)
#' n <- 500
#' d <- data.frame(x = rnorm(n), z = rnorm(n))
#' d$y <- ifelse(rbinom(n, 1, plogis(-1 + 0.8 * d$z)) == 1, 0,
#'               rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
#' fit <- zinbreg(y ~ x | z, data = d)
#' fit
#' summary(fit)
#' coef(fit, model = "zero")
#' AIC(fit)
#' nd <- data.frame(x = c(0, 1), z = c(0, 1))
#' predict(fit, newdata = nd, type = "response")
#' predict(fit, newdata = nd, type = "zero")
#' @name zinbreg-methods
NULL

#' @rdname zinbreg-methods
#' @export
print.zinbreg <- function(x, digits = max(3L, getOption("digits") - 3L),
                          ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Zero-inflated negative binomial regression (EM with IRLS M-steps)\n\n")
  cat("Count model coefficients (negbin with log link):\n")
  print(format(x$coefficients$count, digits = digits), quote = FALSE)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  print(format(x$coefficients$zero, digits = digits), quote = FALSE)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Log-likelihood:", format(x$loglik, digits = digits + 2), "on",
      length(unlist(x$coefficients)) + 1L, "df\n")
  if (!x$converged) cat("Warning: the EM algorithm did not converge\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
summary.zinbreg <- function(object, ...) {
  se <- sqrt(diag(object$vcov_full))
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  # Count table includes Log(theta), as in pscl's summary.
  count <- coef_table(c(object$coefficients$count,
                        `Log(theta)` = log(object$theta)),
                      se[c(seq_len(p), p + q + 1L)])
  zero <- coef_table(object$coefficients$zero, se[p + seq_len(q)])
  structure(list(
    call = object$call,
    coefficients = list(count = count, zero = zero),
    theta = object$theta,
    SE.theta = object$SE.theta,
    loglik = object$loglik,
    df = p + q + 1L,
    aic = stats::AIC(object),
    iter = object$iter,
    converged = object$converged,
    residuals = stats::residuals(object, type = "pearson")
  ), class = "summary.zinbreg")
}

#' @rdname zinbreg-methods
#' @export
print.summary.zinbreg <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("Pearson residuals:\n")
  print(stats::setNames(stats::quantile(x$residuals),
                        c("Min", "1Q", "Median", "3Q", "Max")),
        digits = digits)
  cat("\nCount model coefficients (negbin with log link):\n")
  stats::printCoefmat(x$coefficients$count, digits = digits, ...)
  cat("\nZero-inflation model coefficients (binomial with logit link):\n")
  stats::printCoefmat(x$coefficients$zero, digits = digits, ...)
  cat("\nTheta:", format(x$theta, digits = digits),
      "  Std. Error:", format(x$SE.theta, digits = digits), "\n")
  cat("Log-likelihood:", format(x$loglik, digits = digits + 2), "on",
      x$df, "df   AIC:", format(x$aic, digits = digits + 2), "\n")
  cat("EM iterations:", x$iter,
      if (x$converged) "(converged)" else "(NOT converged)", "\n")
  invisible(x)
}

#' @rdname zinbreg-methods
#' @export
coef.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  switch(model,
         full = stats::setNames(
           unlist(object$coefficients, use.names = FALSE),
           rownames(object$vcov)),
         count = object$coefficients$count,
         zero = object$coefficients$zero)
}

#' @rdname zinbreg-methods
#' @export
vcov.zinbreg <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  p <- length(object$coefficients$count)
  q <- length(object$coefficients$zero)
  idx <- switch(model, full = seq_len(p + q), count = seq_len(p),
                zero = p + seq_len(q))
  object$vcov[idx, idx, drop = FALSE]
}

#' @rdname zinbreg-methods
#' @export
logLik.zinbreg <- function(object, ...) {
  make_logLik(object$loglik, df = length(unlist(object$coefficients)) + 1L,
              nobs = object$nobs)
}

#' @rdname zinbreg-methods
#' @export
nobs.zinbreg <- function(object, ...) object$nobs

#' @rdname zinbreg-methods
#' @export
fitted.zinbreg <- function(object, ...) object$fitted.values

#' @rdname zinbreg-methods
#' @export
residuals.zinbreg <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  r <- object$y - object$fitted.values
  if (type == "pearson") {
    # ZINB variance: E(Y) = (1 - pi) mu and
    # E(Y^2) = (1 - pi)(mu + mu^2 / theta + mu^2), so
    # Var(Y) = (1 - pi) mu (1 + mu (pi + 1 / theta)).
    mu <- object$mu
    pi <- object$pi
    r <- r / sqrt((1 - pi) * mu * (1 + mu * (pi + 1 / object$theta)))
  }
  r
}

#' @rdname zinbreg-methods
#' @export
predict.zinbreg <- function(object, newdata = NULL,
                            type = c("response", "count", "zero"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    mu <- object$mu
    pi <- object$pi
  } else {
    dx <- new_design(object$terms$count, newdata, object$xlevels$count,
                     object$contrasts$count, object$call)
    dz <- new_design(object$terms$zero, newdata, object$xlevels$zero,
                     object$contrasts$zero)
    # mu = exp(x' beta + o),  pi = logit^{-1}(z' gamma)
    mu <- exp(drop(dx$X %*% object$coefficients$count) + dx$offset)
    pi <- stats::plogis(drop(dz$X %*% object$coefficients$zero))
    names(mu) <- names(pi) <- rownames(newdata)
  }
  switch(type,
         response = (1 - pi) * mu,
         count = mu,
         zero = pi)
}
