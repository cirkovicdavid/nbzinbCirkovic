# Zero-inflated negative binomial regression

Fits a zero-inflated negative binomial (ZINB) model. With probability
\\\pi_i\\ observation \\i\\ is a structural zero; otherwise it comes
from an NB2 distribution with mean \\\mu_i\\ and dispersion \\\theta\\.
The count model is \\\log\mu_i = x_i^\top\beta + o_i\\ and the
zero-inflation model is \\\mathrm{logit}\\\pi_i = z_i^\top\gamma\\.

## Usage

``` r
zinbreg(formula, data, offset, control = zinbreg_control())
```

## Arguments

- formula:

  A two-part formula `y ~ x1 + x2 | z1 + z2`: count-model covariates
  before `|`, zero-model covariates after it. Without `|`, the same
  covariates are used in both parts. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term is allowed in
  the count part only.

- data:

  A data frame containing the variables in `formula`.

- offset:

  An optional offset for the count model (e.g. log exposure), evaluated
  in `data`.

- control:

  A list of settings from
  [`zinbreg_control()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/zinbreg_control.md).

## Value

An object of class `"zinbreg"`, a list with components

- coefficients:

  a list with elements `count` (\\\hat\beta\\) and `zero`
  (\\\hat\gamma\\)

- vcov:

  variance-covariance matrix of \\(\hat\beta, \hat\gamma)\\; names are
  prefixed `count_` and `zero_`

- theta, SE.theta, SE.logtheta:

  estimated dispersion parameter and the standard errors of
  \\\hat\theta\\ and \\\log\hat\theta\\

- vcov_full:

  variance-covariance matrix of \\(\hat\beta, \hat\gamma,
  \log\hat\theta)\\

- loglik, loglik_path:

  maximized log-likelihood and its value at every EM iteration

- iter, converged:

  number of EM iterations and convergence flag

- fitted.values:

  \\(1 - \hat\pi_i)\hat\mu_i\\, the fitted mean of \\Y_i\\

- mu, pi, tau:

  fitted count means, structural-zero probabilities, and final E-step
  posterior probabilities

- y, offset, nobs, df.residual:

  response, offset, number of observations, and residual degrees of
  freedom

- call, formula, terms, model, xlevels, contrasts, control:

  information used to print the fit and to predict for new data

## Details

The model is fitted by the EM algorithm. The E-step computes the
posterior probability \\\tau_i\\ that each zero is a structural zero.
The M-step fits (i) a logistic regression of \\\tau_i\\ on \\z_i\\ by
IRLS and (ii) a weighted NB regression with prior weights \\1 - \tau_i\\
by the same IRLS / Newton-Raphson algorithm as
[`nbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/nbreg.md).
Iterations stop when the observed-data log-likelihood changes by less
than `control$epsilon`. No general-purpose numerical maximizer is used.

EM does not produce standard errors. They are obtained by inverting the
negative Hessian of the observed-data log-likelihood with respect to
\\(\beta, \gamma, \log\theta)\\ at the estimates. By default the Hessian
is computed from analytic second derivatives;
`zinbreg_control(hessian = "numeric")` uses
[`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)
instead.

## See also

[`zinbreg_control()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/zinbreg_control.md),
[`nbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/nbreg.md)

## Examples

``` r
# Simulated ZINB data
set.seed(1)
n <- 500
x <- rnorm(n)
z <- rnorm(n)
structural_zero <- rbinom(n, 1, plogis(-1 + 0.8 * z))
y <- ifelse(structural_zero == 1, 0,
            rnbinom(n, size = 2, mu = exp(1 + 0.5 * x)))
fit <- zinbreg(y ~ x | z, data = data.frame(x, z, y))
fit$coefficients
#> $count
#> (Intercept)           x 
#>   1.0953534   0.4391556 
#> 
#> $zero
#> (Intercept)           z 
#>  -1.0497237   0.8567723 
#> 
fit$theta
#> [1] 1.960634

# Publication counts of biochemistry PhD students
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  fit2 <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                  data = bioChemists)
  fit2$coefficients
}
#> $count
#>  (Intercept)     femWomen   marMarried         kid5          phd         ment 
#>  0.404013648 -0.211906301  0.139462144 -0.167623779  0.001962561  0.024392582 
#> 
#> $zero
#> (Intercept)        ment 
#>  -0.8066567  -0.6094417 
#> 
```
