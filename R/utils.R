# Shared helpers for the S3 methods.

# Build the model matrix and offset for new data, using the terms, factor
# levels (xlevels) and contrasts saved at fitting time, so that new data
# are coded exactly like the original data.
#   terms      terms object of the fitted model part (response is dropped)
#   newdata    data frame of new covariate values
#   xlevels    factor levels saved by .getXlevels() when fitting
#   contrasts  contrasts saved from the fitted model matrix
#   call       the original call (to re-evaluate an offset = argument)
new_design <- function(terms, newdata, xlevels, contrasts, call = NULL) {
  tt <- stats::delete.response(terms)
  mf <- stats::model.frame(tt, newdata, na.action = stats::na.pass,
                           xlev = xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  # Offsets: offset() terms in the formula ...
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(X))
  # ... plus an offset = argument, evaluated in the new data.
  if (!is.null(call) && !is.null(call$offset)) {
    off <- off + eval(call$offset, newdata, environment(terms))
  }
  list(X = X, offset = off)
}

# Coefficient table used by the summary() methods: estimate, standard error,
# Wald z = estimate / SE, and two-sided p-value 2 P(Z > |z|).
coef_table <- function(est, se) {
  z <- est / se
  cbind(Estimate = est, `Std. Error` = se, `z value` = z,
        `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
}

# Build a "logLik" object; the df attribute makes AIC() and BIC() work, and
# nobs is used by BIC().
make_logLik <- function(value, df, nobs) {
  structure(value, df = df, nobs = nobs, class = "logLik")
}
