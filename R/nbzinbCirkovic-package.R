#' nbzinbCirkovic: Negative Binomial and Zero-Inflated Negative Binomial
#' Regression
#'
#' Fits NB2 regression ([nbreg()]) by IRLS alternated with Newton-Raphson for
#' \eqn{\log\theta}{log(theta)}, and zero-inflated NB regression
#' ([zinbreg()]) by an EM algorithm whose M-steps are IRLS fits. Standard S3
#' methods (`summary`, `predict`, `residuals`, `logLik`, ...) are provided
#' for both.
#'
#' @keywords internal
#' @importFrom stats coef vcov logLik nobs fitted residuals predict AIC BIC
"_PACKAGE"
