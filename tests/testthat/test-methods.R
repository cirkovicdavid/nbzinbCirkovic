# S3 methods (Section 4.2): compare with the reference implementations.

nb_pair <- function() {
  data("bioChemists", package = "pscl")
  f <- art ~ fem + mar + kid5 + phd + ment
  list(fit = nbreg(f, data = bioChemists),
       ref = MASS::glm.nb(f, data = bioChemists,
                          control = stats::glm.control(epsilon = 1e-10,
                                                       maxit = 100)),
       data = bioChemists)
}

zinb_pair <- function() {
  data("bioChemists", package = "pscl")
  f <- art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  list(fit = zinbreg(f, data = bioChemists),
       ref = pscl::zeroinfl(f, data = bioChemists, dist = "negbin",
                            control = pscl::zeroinfl.control(reltol = 1e-12)),
       data = bioChemists)
}

test_that("nbreg methods agree with glm.nb", {
  skip_if_not_installed("MASS")
  skip_if_not_installed("pscl")
  m <- nb_pair()
  fit <- m$fit
  ref <- m$ref
  expect_equal(coef(fit), coef(ref), tolerance = 1e-6)
  expect_equal(vcov(fit), vcov(ref), tolerance = 1e-6)
  expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
  expect_equal(AIC(fit), AIC(ref), tolerance = 1e-8)
  expect_equal(BIC(fit), BIC(ref), tolerance = 1e-8)
  expect_equal(nobs(fit), nrow(m$data))
  expect_equal(fitted(fit), fitted(ref), tolerance = 1e-6)
  expect_equal(residuals(fit, "pearson"), residuals(ref, "pearson"),
               tolerance = 1e-6)
  expect_equal(residuals(fit, "response"), residuals(ref, "response"),
               tolerance = 1e-6)
  expect_equal(predict(fit, type = "link"), predict(ref, type = "link"),
               tolerance = 1e-6)
  # Wald table matches summary.glm
  s <- summary(fit)
  expect_s3_class(s, "summary.nbreg")
  expect_equal(unname(s$coefficients), unname(coef(summary(ref))),
               tolerance = 1e-5)
})

test_that("nbreg predict() with newdata reproduces in-sample predictions", {
  skip_if_not_installed("pscl")
  m <- nb_pair()
  # A subset in which the factors do not have all their levels must still
  # be coded with the original levels.
  nd <- m$data[m$data$fem == "Women" & m$data$mar == "Married", ][1:5, ]
  expect_equal(predict(m$fit, newdata = nd),
               fitted(m$fit)[rownames(nd)])
  expect_equal(predict(m$fit, newdata = nd, type = "link"),
               log(fitted(m$fit)[rownames(nd)]))
})

test_that("nbreg predict() applies the offset to new data", {
  set.seed(3)
  d <- data.frame(x = rnorm(200), t = runif(200, 1, 4))
  d$y <- rnbinom(200, size = 2, mu = d$t * exp(0.3 + 0.5 * d$x))
  fit_arg <- nbreg(y ~ x, data = d, offset = log(t))
  fit_frm <- nbreg(y ~ x + offset(log(t)), data = d)
  nd <- data.frame(x = c(0, 1), t = c(1, 2))
  expected <- nd$t * exp(coef(fit_arg)[1] + coef(fit_arg)[2] * nd$x)
  expect_equal(unname(predict(fit_arg, newdata = nd)), unname(expected))
  expect_equal(unname(predict(fit_frm, newdata = nd)), unname(expected))
})

test_that("zinbreg methods agree with zeroinfl", {
  skip_if_not_installed("pscl")
  m <- zinb_pair()
  fit <- m$fit
  ref <- m$ref
  expect_equal(names(coef(fit)), names(coef(ref)))
  expect_lt(max(abs(coef(fit) - coef(ref))), 1e-3)
  expect_equal(coef(fit, "count"), coef(ref, "count"), tolerance = 1e-3)
  expect_equal(coef(fit, "zero"), coef(ref, "zero"), tolerance = 1e-3)
  expect_equal(dim(vcov(fit, "zero")), dim(vcov(ref, "zero")))
  expect_equal(attr(logLik(fit), "df"), attr(logLik(ref), "df"))
  expect_lt(abs(AIC(fit) - AIC(ref)), 1e-2)
  expect_lt(abs(BIC(fit) - BIC(ref)), 1e-2)
  expect_equal(nobs(fit), nrow(m$data))  # zeroinfl has no nobs() method
  for (type in c("response", "count", "zero")) {
    expect_lt(max(abs(predict(fit, type = type) -
                        predict(ref, type = type))), 1e-3)
  }
  expect_lt(max(abs(fitted(fit) - fitted(ref))), 1e-3)
  expect_lt(max(abs(residuals(fit, "pearson") -
                      residuals(ref, "pearson"))), 1e-3)
  expect_lt(max(abs(residuals(fit, "response") -
                      residuals(ref, "response"))), 1e-3)
  # summary tables, including the Log(theta) row
  s <- summary(fit)
  sr <- summary(ref)
  expect_equal(rownames(s$coefficients$count),
               rownames(sr$coefficients$count))
  expect_lt(max(abs(s$coefficients$count[, 1:2] -
                      sr$coefficients$count[, 1:2])), 1e-3)
  expect_lt(max(abs(s$coefficients$zero[, 1:2] -
                      sr$coefficients$zero[, 1:2])), 1e-3)
})

test_that("zinbreg predict() with newdata reproduces in-sample predictions", {
  skip_if_not_installed("pscl")
  m <- zinb_pair()
  nd <- m$data[m$data$fem == "Men", ][1:4, ]
  for (type in c("response", "count", "zero")) {
    expect_equal(predict(m$fit, newdata = nd, type = type),
                 predict(m$fit, type = type)[rownames(nd)])
  }
})

test_that("print and summary methods print without error", {
  set.seed(4)
  d <- data.frame(x = rnorm(300), z = rnorm(300))
  d$y <- ifelse(rbinom(300, 1, plogis(-1 + d$z)) == 1, 0,
                rnbinom(300, size = 2, mu = exp(0.8 + 0.4 * d$x)))
  nb <- nbreg(y ~ x, data = d)
  zi <- zinbreg(y ~ x | z, data = d)
  expect_output(print(nb), "Theta")
  expect_output(print(summary(nb)), "z value")
  expect_output(print(zi), "Zero-inflation")
  expect_output(print(summary(zi)), "Log\\(theta\\)")
})
