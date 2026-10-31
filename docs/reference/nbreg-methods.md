# Methods for nbreg fits

Standard methods for objects returned by
[`nbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/nbreg.md).

## Usage

``` r
# S3 method for class 'nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
summary(object, ...)

# S3 method for class 'summary.nbreg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'nbreg'
coef(object, ...)

# S3 method for class 'nbreg'
vcov(object, ...)

# S3 method for class 'nbreg'
logLik(object, ...)

# S3 method for class 'nbreg'
nobs(object, ...)

# S3 method for class 'nbreg'
fitted(object, ...)

# S3 method for class 'nbreg'
residuals(object, type = c("pearson", "response"), ...)

# S3 method for class 'nbreg'
predict(object, newdata = NULL, type = c("response", "link"), ...)
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

  An object of class `"nbreg"` (or `"summary.nbreg"` for the `print`
  method of a summary).

- type:

  For [`residuals()`](https://rdrr.io/r/stats/residuals.html):
  `"pearson"` (default), \\(y_i - \hat\mu_i) / \sqrt{\hat\mu_i +
  \hat\mu_i^2/\hat\theta}\\, or `"response"`, \\y_i - \hat\mu_i\\. For
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"response"`
  (default) for the mean \\\mu\\, or `"link"` for the linear predictor
  \\\eta = \log\mu\\.

- newdata:

  Optional data frame of covariate values for prediction. If omitted,
  predictions for the data used in the fit are returned.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html) the coefficient vector;
[`vcov()`](https://rdrr.io/r/stats/vcov.html) its covariance matrix;
[`logLik()`](https://rdrr.io/r/stats/logLik.html) a `"logLik"` object
with `df` = number of coefficients + 1 (for \\\theta\\);
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) numeric vectors;
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the number of
observations; [`summary()`](https://rdrr.io/r/base/summary.html) an
object of class `"summary.nbreg"` with the coefficient table.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200))
d$y <- rnbinom(200, size = 2, mu = exp(0.5 + 0.3 * d$x))
fit <- nbreg(y ~ x, data = d)
fit
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Negative binomial regression (IRLS + Newton-Raphson for theta)
#> 
#> Coefficients:
#> (Intercept)           x 
#>      0.5139      0.3306 
#> 
#> Theta: 1.804   Log-likelihood: -352.486 on 3 df
summary(fit)
#> 
#> Call:
#> nbreg(formula = y ~ x, data = d)
#> 
#> Pearson residuals:
#>     Min      1Q  Median      3Q     Max 
#> -1.1012 -0.8336 -0.2219  0.3826  4.7518 
#> 
#> Coefficients (negative binomial with log link):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  0.51389    0.07705   6.670 2.56e-11 ***
#> x            0.33060    0.08135   4.064 4.82e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Theta: 1.804   Std. Error: 0.4004 
#> Log-likelihood: -352.486 on 3 df   AIC: 710.971 
#> Iterations: 2 (converged) 
AIC(fit)
#> [1] 710.9714
head(residuals(fit, type = "pearson"))
#>          1          2          3          4          5          6 
#> 1.71091893 1.18427336 0.49794105 2.65622086 0.06975723 0.49185960 
predict(fit, newdata = data.frame(x = c(-1, 0, 1)))
#>        1        2        3 
#> 1.201162 1.671786 2.326805 
```
