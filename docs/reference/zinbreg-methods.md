# Methods for zinbreg fits

Standard methods for objects returned by
[`zinbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/zinbreg.md).

## Usage

``` r
# S3 method for class 'zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
summary(object, ...)

# S3 method for class 'summary.zinbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'zinbreg'
coef(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
vcov(object, model = c("full", "count", "zero"), ...)

# S3 method for class 'zinbreg'
logLik(object, ...)

# S3 method for class 'zinbreg'
nobs(object, ...)

# S3 method for class 'zinbreg'
fitted(object, ...)

# S3 method for class 'zinbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'zinbreg'
predict(object, newdata = NULL, type = c("response", "count", "zero"), ...)
```

## Arguments

- digits:

  Number of significant digits to print.

- ...:

  Further arguments (currently ignored, except by
  [`print()`](https://rdrr.io/r/base/print.html) methods which pass them
  to
  [`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html)).

- object, x:

  An object of class `"zinbreg"` (or `"summary.zinbreg"` for the `print`
  method of a summary).

- model:

  For [`coef()`](https://rdrr.io/r/stats/coef.html) and
  [`vcov()`](https://rdrr.io/r/stats/vcov.html): `"full"` (default;
  count and zero coefficients, named `count_*` and `zero_*`), `"count"`,
  or `"zero"`.

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"pearson"` (default) or `"response"`. The response residual is
  \\y_i - (1 - \hat\pi_i)\hat\mu_i\\; the Pearson residual divides it by
  \\\sqrt{(1-\pi_i)\mu_i(1 + \mu_i(\pi_i + 1/\theta))}\\, the standard
  deviation of a ZINB variable. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (default) for the mean \\(1 - \pi)\mu\\, `"count"` for the count-model
  mean \\\mu\\, or `"zero"` for the structural-zero probability \\\pi\\.

- newdata:

  Optional data frame of covariate values for prediction.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) the coefficients and their
covariance matrix; [`logLik()`](https://rdrr.io/r/stats/logLik.html) a
`"logLik"` object with `df` = number of count and zero coefficients + 1
(for \\\theta\\);
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) numeric vectors;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
observations; [`summary()`](https://rdrr.io/r/base/summary.html) an
object of class `"summary.zinbreg"`.

## Examples

``` r
set.seed(1)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$y <- ifelse(rbinom(n, 1, plogis(-1 + 0.8 * d$z)) == 1, 0,
              rnbinom(n, size = 2, mu = exp(1 + 0.5 * d$x)))
fit <- zinbreg(y ~ x | z, data = d)
fit
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Zero-inflated negative binomial regression (EM with IRLS M-steps)
#> 
#> Count model coefficients (negbin with log link):
#> (Intercept)           x 
#>      1.0954      0.4392 
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#> (Intercept)           z 
#>     -1.0497      0.8568 
#> 
#> Theta: 1.961   Log-likelihood: -956.838 on 5 df
summary(fit)
#> 
#> Call:
#> zinbreg(formula = y ~ x | z, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.0637 -0.7140 -0.3876  0.5721  5.7168 
#> 
#> Count model coefficients (negbin with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.09535    0.06121  17.895  < 2e-16 ***
#> x            0.43916    0.05255   8.357  < 2e-16 ***
#> Log(theta)   0.67327    0.19050   3.534 0.000409 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Zero-inflation model coefficients (binomial with logit link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -1.0497     0.1972  -5.324 1.01e-07 ***
#> z             0.8568     0.1570   5.458 4.81e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta: 1.961   Std. Error: 0.3735 
#> Log-likelihood: -956.838 on 5 df   AIC: 1923.68 
#> EM iterations: 34 (converged) 
coef(fit, model = "zero")
#> (Intercept)           z 
#>  -1.0497237   0.8567723 
AIC(fit)
#> [1] 1923.676
nd <- data.frame(x = c(0, 1), z = c(0, 1))
predict(fit, newdata = nd, type = "response")
#>        1        2 
#> 2.214935 2.542610 
predict(fit, newdata = nd, type = "zero")
#>         1         2 
#> 0.2592782 0.4519112 
```
