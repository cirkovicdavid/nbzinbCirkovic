# Compare nbreg() with MASS::glm.nb() using the tolerances of Table 4:
# coefficients and SEs <= 1e-5, theta and log-likelihood <= 1e-4.
expect_matches_glm_nb <- function(fit, ref) {
  expect_true(fit$converged)
  expect_equal(names(fit$coefficients), names(coef(ref)))
  expect_lt(max(abs(fit$coefficients - coef(ref))), 1e-5)
  expect_lt(max(abs(sqrt(diag(fit$vcov)) - sqrt(diag(vcov(ref))))), 1e-5)
  expect_lt(abs(fit$theta - ref$theta), 1e-4)
  expect_lt(abs(fit$SE.theta - ref$SE.theta), 1e-4)
  expect_lt(abs(fit$loglik - as.numeric(logLik(ref))), 1e-4)
}

ref_control <- function() {
  stats::glm.control(epsilon = 1e-10, maxit = 100)
}

test_that("nbreg matches glm.nb on bioChemists", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  data("bioChemists", package = "pscl")
  f <- art ~ fem + mar + kid5 + phd + ment
  fit <- nbreg(f, data = bioChemists)
  ref <- MASS::glm.nb(f, data = bioChemists, control = ref_control())
  expect_matches_glm_nb(fit, ref)
  expect_s3_class(fit, "nbreg")
  expect_equal(fit$nobs, nrow(bioChemists))
})

test_that("nbreg matches glm.nb on NMES1988", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("AER")
  data("NMES1988", package = "AER")
  f <- visits ~ hospital + health + chronic + gender + school + insurance
  fit <- nbreg(f, data = NMES1988)
  ref <- MASS::glm.nb(f, data = NMES1988, control = ref_control())
  expect_matches_glm_nb(fit, ref)
})

# Simulated data with an exposure variable, used by the next tests.
sim_nb <- function(n = 400, seed = 651) {
  set.seed(seed)
  x <- rnorm(n)
  g <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  t <- runif(n, 0.5, 3)
  mu <- t * exp(0.2 + 0.4 * x + 0.3 * (g == "b") - 0.5 * (g == "c"))
  data.frame(y = rnbinom(n, size = 1.5, mu = mu), x = x, g = g, t = t)
}

test_that("offset argument and offset() in the formula agree with glm.nb", {
  skip_if_not_installed("MASS")
  d <- sim_nb()
  fit_arg <- nbreg(y ~ x + g, data = d, offset = log(t))
  fit_frm <- nbreg(y ~ x + g + offset(log(t)), data = d)
  ref <- MASS::glm.nb(y ~ x + g + offset(log(t)), data = d,
                      control = ref_control())
  expect_matches_glm_nb(fit_arg, ref)
  expect_equal(fit_arg$coefficients, fit_frm$coefficients)
  expect_equal(fit_arg$offset, log(d$t))
})

test_that("single-covariate and intercept-only models work", {
  skip_if_not_installed("MASS")
  d <- sim_nb()
  fit1 <- nbreg(y ~ x, data = d)
  ref1 <- MASS::glm.nb(y ~ x, data = d, control = ref_control())
  expect_matches_glm_nb(fit1, ref1)

  # Intercept only: the MLE of the intercept is log(mean(y)).
  fit0 <- nbreg(y ~ 1, data = d)
  expect_equal(unname(fit0$coefficients), log(mean(d$y)), tolerance = 1e-8)
})

test_that("nbreg gives informative errors for bad input", {
  d <- sim_nb(n = 50)
  d_neg <- d
  d_neg$y[1] <- -1
  expect_error(nbreg(y ~ x, data = d_neg), "non-negative integer")

  d_frac <- d
  d_frac$y[1] <- 1.5
  expect_error(nbreg(y ~ x, data = d_frac), "non-negative integer")

  d$x2 <- 2 * d$x  # perfectly collinear with x
  expect_error(nbreg(y ~ x + x2, data = d), "rank deficient")

  expect_error(nbreg(y ~ x, data = d, control = nbreg_control(epsilon = -1)),
               "epsilon")
  expect_error(nbreg(y ~ x, data = d, control = nbreg_control(maxit = 0)),
               "maxit")
})

test_that("nbreg warns when it does not converge", {
  d <- sim_nb()
  expect_warning(nbreg(y ~ x, data = d, control = nbreg_control(maxit = 1)),
                 "did not converge")
})

test_that("nbreg warns when theta is very large (Poisson data)", {
  set.seed(1)
  x <- rnorm(500)
  d <- data.frame(x = x, y = rpois(500, exp(1 + 0.3 * x)))
  expect_warning(fit <- nbreg(y ~ x, data = d), "theta estimate is very large")
  expect_gt(fit$theta, 1e4)
})

test_that("trace = TRUE prints the log-likelihood each iteration", {
  d <- sim_nb(n = 100)
  expect_message(nbreg(y ~ x, data = d, control = nbreg_control(trace = TRUE)),
                 "logLik")
})
