# The analytic score and Hessian (R/hessian.R) must agree with numerical
# derivatives of our own log-likelihood, at the MLE and away from it.
deriv_check <- function(fit, shift_sd = 0.1, seed = 2) {
  X <- model.matrix(fit$terms$count, fit$model)
  Z <- model.matrix(fit$terms$zero, fit$model)
  y <- fit$y
  o <- fit$offset
  par_mle <- c(fit$coefficients$count, fit$coefficients$zero, log(fit$theta))
  set.seed(seed)
  par_off <- par_mle + rnorm(length(par_mle), sd = shift_sd)
  for (par in list(par_mle, par_off)) {
    H_a <- zinb_hessian(par, X, Z, y, o)
    H_n <- numDeriv::hessian(zinb_loglik_par, par, X = X, Z = Z, y = y,
                             offset = o)
    expect_lt(max(abs(H_a - H_n)) / max(abs(H_n)), 1e-6)
    expect_true(isSymmetric(unname(H_a)))
    g_a <- zinb_score(par, X, Z, y, o)
    g_n <- numDeriv::grad(zinb_loglik_par, par, X = X, Z = Z, y = y,
                          offset = o)
    expect_lt(max(abs(g_a - g_n)), 1e-3 * max(1, max(abs(g_n))))
  }
}

test_that("analytic derivatives match numDeriv on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  deriv_check(fit)
})

test_that("analytic derivatives match numDeriv with an offset", {
  set.seed(7)
  n <- 400
  d <- data.frame(x = rnorm(n), z = rnorm(n), t = runif(n, 0.5, 3))
  d$y <- ifelse(rbinom(n, 1, plogis(-0.5 + d$z)) == 1, 0,
                rnbinom(n, size = 1.5, mu = d$t * exp(0.5 + 0.3 * d$x)))
  fit <- zinbreg(y ~ x | z, data = d, offset = log(t))
  deriv_check(fit)
})

test_that("analytic and numeric Hessians give the same standard errors", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  f <- art ~ fem + kid5 + ment | ment
  fit_a <- zinbreg(f, data = bioChemists)
  fit_n <- zinbreg(f, data = bioChemists,
                   control = zinbreg_control(hessian = "numeric"))
  expect_equal(fit_a$control$hessian, "analytic")
  expect_lt(max(abs(sqrt(diag(fit_a$vcov_full)) -
                      sqrt(diag(fit_n$vcov_full)))), 1e-5)
})

test_that("the score is close to zero at the EM estimates", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists,
                 control = zinbreg_control(epsilon = 1e-12))
  X <- model.matrix(fit$terms$count, fit$model)
  Z <- model.matrix(fit$terms$zero, fit$model)
  par <- c(fit$coefficients$count, fit$coefficients$zero, log(fit$theta))
  expect_lt(max(abs(zinb_score(par, X, Z, fit$y, fit$offset))), 1e-3)
})
