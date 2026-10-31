# Control settings for nbreg()

Control settings for nbreg()

## Usage

``` r
nbreg_control(
  epsilon = 1e-08,
  maxit = 100,
  irls_maxit = 100,
  theta_maxit = 100,
  theta_warn = 10000,
  trace = FALSE
)
```

## Arguments

- epsilon:

  Convergence tolerance. The outer loop stops when the absolute change
  in the log-likelihood is below `epsilon`; the inner IRLS and
  Newton-Raphson loops stop when their parameter updates are below
  `epsilon`.

- maxit:

  Maximum number of outer iterations (each one is a Newton-Raphson
  update of theta followed by an IRLS fit of beta).

- irls_maxit:

  Maximum number of IRLS iterations per outer iteration.

- theta_maxit:

  Maximum number of Newton-Raphson iterations per outer iteration.

- theta_warn:

  Warn if the estimated theta exceeds this value (the Poisson limit).

- trace:

  Logical; if `TRUE`, print the log-likelihood and theta at each outer
  iteration.

## Value

A list of control settings.

## Examples

``` r
nbreg_control(epsilon = 1e-10, trace = TRUE)
#> $epsilon
#> [1] 1e-10
#> 
#> $maxit
#> [1] 100
#> 
#> $irls_maxit
#> [1] 100
#> 
#> $theta_maxit
#> [1] 100
#> 
#> $theta_warn
#> [1] 10000
#> 
#> $trace
#> [1] TRUE
#> 
```
