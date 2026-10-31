# Control settings for zinbreg()

Control settings for zinbreg()

## Usage

``` r
zinbreg_control(
  epsilon = 1e-08,
  maxit = 5000,
  irls_maxit = 100,
  theta_maxit = 100,
  nb_maxit = 100,
  theta_warn = 10000,
  hessian = c("analytic", "numeric"),
  trace = FALSE
)
```

## Arguments

- epsilon:

  Convergence tolerance. EM stops when the absolute change in the
  observed-data log-likelihood is below `epsilon`. The same tolerance is
  used inside the M-steps.

- maxit:

  Maximum number of EM iterations.

- irls_maxit:

  Maximum number of IRLS iterations in each M-step fit.

- theta_maxit:

  Maximum number of Newton-Raphson iterations for theta in each count
  M-step.

- nb_maxit:

  Maximum number of outer (theta / beta) iterations in each count
  M-step.

- theta_warn:

  Warn if the estimated theta exceeds this value.

- hessian:

  How to compute the Hessian for the standard errors: `"analytic"`
  (default, exact second derivatives) or `"numeric"` (finite differences
  with
  [`numDeriv::hessian()`](https://rdrr.io/pkg/numDeriv/man/hessian.html)).

- trace:

  Logical; if `TRUE`, print the observed-data log-likelihood and its
  change at each EM iteration.

## Value

A list of control settings.

## Examples

``` r
zinbreg_control(epsilon = 1e-10, trace = TRUE)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
#> [1] 5000
#> 
#> $irls_maxit
#> [1] 100
#> 
#> $theta_maxit
#> [1] 100
#> 
#> $nb_maxit
#> [1] 100
#> 
#> $theta_warn
#> [1] 10000
#> 
#> $hessian
#> [1] "analytic"
#> 
#> $trace
#> [1] TRUE
#> 
zinbreg_control(hessian = "numeric")
#> $epsilon
#> [1] 1e-08
#> 
#> $maxit
#> [1] 5000
#> 
#> $irls_maxit
#> [1] 100
#> 
#> $theta_maxit
#> [1] 100
#> 
#> $nb_maxit
#> [1] 100
#> 
#> $theta_warn
#> [1] 10000
#> 
#> $hessian
#> [1] "numeric"
#> 
#> $trace
#> [1] FALSE
#> 
```
