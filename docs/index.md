# nbzinbCirkovic

Negative binomial (NB) and zero-inflated negative binomial (ZINB)
regression, with the estimation implemented from scratch:

- [`nbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/nbreg.md)
  fits NB2 regression by **iteratively reweighted least squares (IRLS)**
  for the coefficients, alternated with **Newton–Raphson** updates of
  log(θ).
- [`zinbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/zinbreg.md)
  fits ZINB regression by an **EM algorithm** whose M-steps are IRLS
  fits (a weighted logistic regression for the zero part and a weighted
  NB regression for the count part). Standard errors come from an
  **analytic Hessian** of the observed-data log-likelihood.

No general-purpose numerical maximizer is used. The results agree with
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) and
[`pscl::zeroinfl()`](https://rdrr.io/pkg/pscl/man/zeroinfl.html); see
the benchmark in the report.

STA 651 final project (Fall 2026), David Cirkovic.

- **Website:** <https://cirkovicdavid.github.io/nbzinbCirkovic/>
  (function reference and the project report)
- **Report (PDF):**
  [`report/report.pdf`](https://cirkovicdavid.github.io/nbzinbCirkovic/report/report.pdf)

## Installation

``` r

# install.packages("devtools")
devtools::install_github("cirkovicdavid/nbzinbCirkovic")
```

Since devtools 2.5.0, `install_github()` needs the `remotes` package; if
the command above reports that “remotes” is required, run
`install.packages("remotes")` first, or use either of

``` r

remotes::install_github("cirkovicdavid/nbzinbCirkovic")
pak::pak("cirkovicdavid/nbzinbCirkovic")
```

## Usage

``` r

library(nbzinbCirkovic)
data("bioChemists", package = "pscl")

# Negative binomial regression
fit_nb <- nbreg(art ~ fem + mar + kid5 + phd + ment, data = bioChemists)
summary(fit_nb)
#> Coefficients (negative binomial with log link):
#>              Estimate Std. Error z value Pr(>|z|)
#> (Intercept)  0.256144   0.137348   1.865 0.062191 .
#> femWomen    -0.216418   0.072636  -2.979 0.002887 **
#> marMarried   0.150489   0.082097   1.833 0.066791 .
#> kid5        -0.176415   0.052813  -3.340 0.000837 ***
#> phd          0.015271   0.035873   0.426 0.670326
#> ment         0.029082   0.003214   9.048  < 2e-16 ***
#>
#> Theta: 2.264   Std. Error: 0.2712
#> Log-likelihood: -1560.96 on 7 df   AIC: 3135.92

# Zero-inflated NB: count covariates | zero-inflation covariates
fit_zinb <- zinbreg(art ~ fem + mar + kid5 + phd + ment | ment,
                    data = bioChemists)
summary(fit_zinb)

# Standard methods work for both classes
AIC(fit_nb, fit_zinb)
predict(fit_zinb, newdata = bioChemists[1:3, ], type = "zero")  # P(structural zero)
residuals(fit_nb, type = "pearson")
```

Use `nbreg_control(trace = TRUE)` or `zinbreg_control(trace = TRUE)` to
print the log-likelihood at every iteration; for
[`zinbreg()`](https://cirkovicdavid.github.io/nbzinbCirkovic/reference/zinbreg.md)
it never decreases, as the EM algorithm guarantees.

## R CMD check

[`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
finishes with 0 errors, 0 warnings, and 0 notes.
