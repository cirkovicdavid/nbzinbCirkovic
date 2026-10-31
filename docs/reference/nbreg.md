# Negative binomial regression

Fits a negative binomial (NB2) regression model with log link, \\\log
\mu_i = x_i^\top \beta + o_i\\ and \\\mathrm{Var}(Y_i) = \mu_i +
\mu_i^2/\theta\\. The coefficients \\\beta\\ are estimated by
iteratively reweighted least squares (IRLS) with \\\theta\\ held fixed,
alternated with Newton-Raphson updates of \\\log\theta\\ with \\\beta\\
held fixed, until the log-likelihood stops changing. No general-purpose
numerical maximizer is used.

## Usage

``` r
nbreg(formula, data, offset, control = nbreg_control())
```

## Arguments

- formula:

  A model formula, `response ~ covariates`. Offsets may be given in the
  formula with [`offset()`](https://rdrr.io/r/stats/offset.html).

- data:

  A data frame containing the variables in `formula`.

- offset:

  An optional offset vector (e.g. log exposure), evaluated in `data`.
  Added to any [`offset()`](https://rdrr.io/r/stats/offset.html) terms
  in the formula.

- control:

  A list of settings from
  [`nbreg_control()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/nbreg_control.md).

## Value

An object of class `"nbreg"`, a list with components

- coefficients:

  estimated regression coefficients \\\hat\beta\\

- vcov:

  variance-covariance matrix of \\\hat\beta\\

- theta, SE.theta:

  estimated dispersion parameter and its standard error

- loglik:

  maximized log-likelihood

- iter, converged:

  number of outer iterations and whether the log-likelihood converged

- fitted.values, linear.predictors:

  \\\hat\mu_i\\ and \\\hat\eta_i\\

- y, offset, nobs, df.residual:

  response, offset, number of observations, and residual degrees of
  freedom

- call, terms, model, xlevels, contrasts, control:

  information used to print the fit and to predict for new data

## Details

The algorithm starts from a Poisson fit (IRLS with \\\theta = \infty\\)
and the moment estimate \\\theta_0 = n / \sum_i (y_i/\mu_i - 1)^2\\.
Standard errors of \\\hat\beta\\ come from \\(X^\top \hat W X)^{-1}\\,
treating \\\theta\\ as known, and the standard error of \\\hat\theta\\
is \\\\-H(\hat\theta)\\^{-1/2}\\, where \\H\\ is the second derivative
of the log-likelihood in \\\theta\\.

## See also

[`nbreg_control()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/nbreg_control.md),
[`zinbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/zinbreg.md)

## Examples

``` r
# Simulated data with known beta = (0.5, 0.3) and theta = 2
set.seed(1)
n <- 300
x <- rnorm(n)
y <- rnbinom(n, size = 2, mu = exp(0.5 + 0.3 * x))
fit <- nbreg(y ~ x, data = data.frame(x, y))
fit$coefficients
#> (Intercept)           x 
#>   0.4154160   0.3302077 
fit$theta
#> [1] 1.70792

# Publication counts of biochemistry PhD students
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  fit2 <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
  cbind(estimate = fit2$coefficients, se = sqrt(diag(fit2$vcov)))
}
#>                estimate          se
#> (Intercept)  0.25614402 0.137347701
#> femWomen    -0.21641842 0.072636088
#> marMarried   0.15048945 0.082096542
#> kid5        -0.17641524 0.052813297
#> phd          0.01527116 0.035872991
#> ment         0.02908234 0.003214319
```
