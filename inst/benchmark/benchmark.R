# Benchmark: nbzinbCirkovic vs. the reference implementations
# (project description, Section 5).
#
# For each data set (pscl::bioChemists, AER::NMES1988) and model (NB, ZINB)
# this script records our estimates and standard errors side by side with
# MASS::glm.nb() / pscl::zeroinfl(), their maximum absolute differences
# (Table 4), and computing times. Results are saved to
# inst/benchmark/results/benchmark.rds and used by the report.
#
# Run from the package root:  Rscript inst/benchmark/benchmark.R

library(nbzinbCirkovic)
library(MASS)
library(pscl)
data("bioChemists", package = "pscl")
data("NMES1988", package = "AER")

specs <- list(
  bioChemists = list(
    data = bioChemists,
    nb = art ~ fem + mar + kid5 + phd + ment,
    zinb = art ~ fem + mar + kid5 + phd + ment | fem + mar + kid5 + phd + ment
  ),
  NMES1988 = list(
    data = NMES1988,
    nb = visits ~ hospital + health + chronic + gender + school + insurance,
    zinb = visits ~ hospital + health + chronic + gender + school + insurance |
      hospital + chronic + gender + school + insurance
  )
)

# Reference fits with tight tolerances, as recommended in Section 5.
ref_nb <- function(f, d) {
  glm.nb(f, data = d, control = glm.control(epsilon = 1e-10, maxit = 100))
}
ref_zinb <- function(f, d) {
  zeroinfl(f, data = d, dist = "negbin",
           control = zeroinfl.control(reltol = 1e-12))
}

# Median elapsed time (seconds) over `reps` runs.
time_it <- function(expr_fun, reps = 10) {
  median(replicate(reps, system.time(expr_fun())[["elapsed"]]))
}

results <- list()
for (ds in names(specs)) {
  s <- specs[[ds]]

  # ---- NB ------------------------------------------------------------------
  ours <- nbreg(s$nb, data = s$data)
  ref <- ref_nb(s$nb, s$data)
  nb_tab <- data.frame(
    term = names(coef(ours)),
    est_ours = coef(ours), est_ref = coef(ref),
    se_ours = sqrt(diag(vcov(ours))), se_ref = sqrt(diag(vcov(ref))),
    row.names = NULL
  )
  nb_tab <- rbind(nb_tab, data.frame(
    term = "theta", est_ours = ours$theta, est_ref = ref$theta,
    se_ours = ours$SE.theta, se_ref = ref$SE.theta))
  nb_diff <- c(
    coefficients = max(abs(coef(ours) - coef(ref))),
    std_errors = max(abs(sqrt(diag(vcov(ours))) - sqrt(diag(vcov(ref))))),
    theta = abs(ours$theta - ref$theta),
    loglik = abs(as.numeric(logLik(ours)) - as.numeric(logLik(ref))))
  nb_time <- c(ours = time_it(function() nbreg(s$nb, data = s$data)),
               ref = time_it(function() ref_nb(s$nb, s$data)))

  # ---- ZINB ----------------------------------------------------------------
  ours_z <- zinbreg(s$zinb, data = s$data)
  ref_z <- ref_zinb(s$zinb, s$data)
  zinb_tab <- data.frame(
    term = names(coef(ours_z)),
    est_ours = coef(ours_z), est_ref = coef(ref_z),
    se_ours = sqrt(diag(vcov(ours_z))), se_ref = sqrt(diag(vcov(ref_z))),
    row.names = NULL
  )
  zinb_tab <- rbind(zinb_tab, data.frame(
    term = "Log(theta)", est_ours = log(ours_z$theta),
    est_ref = log(ref_z$theta), se_ours = ours_z$SE.logtheta,
    se_ref = ref_z$SE.logtheta))
  zinb_diff <- c(
    coefficients = max(abs(coef(ours_z) - coef(ref_z))),
    std_errors = max(abs(sqrt(diag(vcov(ours_z))) -
                           sqrt(diag(vcov(ref_z))))),
    log_theta = abs(log(ours_z$theta) - log(ref_z$theta)),
    loglik = abs(as.numeric(logLik(ours_z)) - as.numeric(logLik(ref_z))))
  zinb_time <- c(ours = time_it(function() zinbreg(s$zinb, data = s$data),
                                reps = 5),
                 ref = time_it(function() ref_zinb(s$zinb, s$data), reps = 5))

  results[[ds]] <- list(
    nb = list(table = nb_tab, max_diff = nb_diff, time = nb_time,
              loglik = c(ours = ours$loglik, ref = as.numeric(logLik(ref))),
              iter = ours$iter),
    zinb = list(table = zinb_tab, max_diff = zinb_diff, time = zinb_time,
                loglik = c(ours = ours_z$loglik,
                           ref = as.numeric(logLik(ref_z))),
                iter = ours_z$iter)
  )
  cat("\n==", ds, "==\nNB max |diff|:\n")
  print(signif(nb_diff, 3))
  cat("ZINB max |diff|:\n")
  print(signif(zinb_diff, 3))
  cat("Timing (s): NB", round(nb_time, 3), " ZINB", round(zinb_time, 3), "\n")
}

results$session <- list(date = Sys.Date(), R = R.version.string,
                        package = as.character(packageVersion("nbzinbCirkovic")))
saveRDS(results, file.path("inst", "benchmark", "results", "benchmark.rds"))
