# Compare zinbreg() with pscl::zeroinfl(dist = "negbin") using the
# tolerances of Table 4: coefficients, SEs, log(theta) and log-likelihood
# all within 1e-3.
expect_matches_zeroinfl <- function(fit, ref) {
  expect_true(fit$converged)
  ours <- unlist(fit$coefficients, use.names = FALSE)
  expect_equal(rownames(fit$vcov), names(coef(ref)))
  expect_lt(max(abs(ours - coef(ref))), 1e-3)
  expect_lt(max(abs(sqrt(diag(fit$vcov)) - sqrt(diag(vcov(ref))))), 1e-3)
  expect_lt(abs(log(fit$theta) - log(ref$theta)), 1e-3)
  expect_lt(abs(fit$SE.logtheta - ref$SE.logtheta), 1e-3)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-3)
}

fit_ref <- function(formula, data) {
  pscl::zeroinfl(formula, data = data, dist = "negbin",
                 control = pscl::zeroinfl.control(reltol = 1e-12))
}

test_that("zinbreg matches zeroinfl on bioChemists", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  fit <- zinbreg(f, data = bioChemists)
  expect_matches_zeroinfl(fit, fit_ref(f, bioChemists))
  expect_s3_class(fit, "zinbreg")
})

test_that("zinbreg matches zeroinfl on NMES1988", {
  skip_if_not_installed("pscl")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER")
  f <- visits ~ hospital + health + chronic + gender + school + insurance |
    hospital + chronic + gender + school + insurance
  fit <- zinbreg(f, data = NMES1988)
  expect_matches_zeroinfl(fit, fit_ref(f, NMES1988))
})

test_that("the EM log-likelihood never decreases", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit <- zinbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  expect_true(all(diff(fit$loglik_path) >= -1e-10))
  expect_equal(fit$loglik, fit$loglik_path[length(fit$loglik_path)])
})

test_that("a formula without | uses the same covariates in both parts", {
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  fit1 <- zinbreg(art ~ fem + ment, data = bioChemists)
  fit2 <- zinbreg(art ~ fem + ment | fem + ment, data = bioChemists)
  expect_equal(fit1$coefficients, fit2$coefficients)
  expect_equal(names(fit1$coefficients$zero), c("(Intercept)", "femWomen",
                                               "ment"))
})

# Simulated ZINB data with an exposure variable.
sim_zinb <- function(n = 600, seed = 651) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  t <- runif(n, 0.5, 3)
  structural <- rbinom(n, 1, plogis(-0.8 + 0.7 * z))
  y <- ifelse(structural == 1, 0,
              rnbinom(n, size = 2, mu = t * exp(0.6 + 0.4 * x)))
  data.frame(y = y, x = x, z = z, t = t)
}

test_that("offset and single-covariate models agree with zeroinfl", {
  skip_if_not_installed("pscl")
  d <- sim_zinb()
  fit <- zinbreg(y ~ x | z, data = d, offset = log(t))
  fit_frm <- zinbreg(y ~ x + offset(log(t)) | z, data = d)
  ref <- fit_ref(y ~ x + offset(log(t)) | z, d)
  expect_matches_zeroinfl(fit, ref)
  expect_equal(fit$coefficients, fit_frm$coefficients)
})

test_that("fitted values are (1 - pi) * mu", {
  d <- sim_zinb(n = 300)
  fit <- zinbreg(y ~ x | z, data = d)
  expect_equal(unname(fit$fitted.values), unname((1 - fit$pi) * fit$mu))
  expect_true(all(fit$tau[d$y > 0] == 0))
})

test_that("zinbreg gives informative errors for bad input", {
  d <- sim_zinb(n = 100)
  d_pos <- d
  d_pos$y <- d_pos$y + 1
  expect_error(zinbreg(y ~ x | z, data = d_pos), "no zeros")

  d_neg <- d
  d_neg$y[1] <- -2
  expect_error(zinbreg(y ~ x | z, data = d_neg), "non-negative integer")

  expect_error(zinbreg(y ~ x | z + offset(log(t)), data = d),
               "count part")
  expect_error(zinbreg(y ~ x | z | x, data = d), "at most two parts")
  expect_error(zinbreg_control(maxit = 0), "maxit")
})

test_that("zinbreg warns when EM does not converge", {
  d <- sim_zinb(n = 300)
  expect_warning(zinbreg(y ~ x | z, data = d,
                         control = zinbreg_control(maxit = 2)),
                 "did not converge")
})

test_that("trace = TRUE prints the log-likelihood at each EM iteration", {
  d <- sim_zinb(n = 200)
  expect_message(zinbreg(y ~ x | z, data = d,
                         control = zinbreg_control(trace = TRUE)),
                 "logLik")
})
