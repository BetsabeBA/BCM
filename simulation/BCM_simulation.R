
# simulation/BCM_simulation.R

rm(list = ls())

suppressPackageStartupMessages({
  library(betareg)
  library(MASS)
  library(optimx)
  library(numDeriv)
})

## Utility: clamp values into (0,1)
.clamp01 <- function(z, eps = .Machine$double.eps) {
  pmin(pmax(z, eps), 1 - eps)
}

## load core functions
source(file.path("R", "calibration_core.R"))


###############################################################
## PARAMETERS (edit these if you want other scenarios)
###############################################################

set.seed(123)

theta_linear <- c(1.3, -1.5, x0 = 1.25,
                  phi = 144, phi0 = 144)

theta   <- theta_linear
x0_true <- theta[3]
link    <- "logit"
n       <- 5
k       <- 3
trials  <- 500
seedNum <- 1   # seedNum= 1 for k=3, and seedNum=2 for k=10

set.seed(10000 + seedNum)

## ---- Output directory (matches BCM_summary.R expectation) ----
out_dir <- file.path("simulation", sprintf("n=%d-%s", n, link))
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}


###############################################################
### STORAGE FOR ALL SIX MODELS
###############################################################

### Linear
linear.coeff    <- rep(NA, trials)
linear.stderr   <- rep(NA, trials)
linear.Criteria <- matrix(NA, nrow = trials, ncol = 3)
linear.coverage <- rep(NA, trials)

### Quadratic
quadratic.coeff    <- rep(NA, trials)
quadratic.stderr   <- rep(NA, trials)
quadratic.Criteria <- matrix(NA, nrow = trials, ncol = 3)
quadratic.coverage <- rep(NA, trials)

### BCM (two-stage beta calibration)
bcm.coeff      <- rep(NA, trials)
bcm.stderr     <- rep(NA, trials)
bcm.Criteria   <- matrix(NA, nrow = trials, ncol = 3)
bcm.coverage   <- rep(NA, trials)
bcm.converge   <- rep(NA, trials)

### WLS
wls.coeff      <- rep(NA, trials)
wls.stderr     <- rep(NA, trials)
wls.Criteria   <- matrix(NA, nrow = trials, ncol = 3)
wls.coverage   <- rep(NA, trials)

### One-stage Beta regression
beta1.coeff    <- rep(NA, trials)
beta1.stderr   <- rep(NA, trials)
beta1.Criteria <- matrix(NA, nrow = trials, ncol = 3)
beta1.coverage <- rep(NA, trials)

### Beta GAM
gam.coeff      <- rep(NA, trials)
gam.stderr     <- rep(NA, trials)
gam.Criteria   <- matrix(NA, nrow = trials, ncol = 3)
gam.coverage   <- rep(NA, trials)

### Neural network
ann.coeff      <- rep(NA, trials)
ann.stderr     <- rep(NA, trials)
ann.Criteria   <- matrix(NA, nrow = trials, ncol = 3)
ann.coverage   <- rep(NA, trials)

## --------------------------------------------
## Convergence tracker (add near your other inits)
## --------------------------------------------
add_conv <- function(df, model, trial, converged, code = NA_integer_, message = NA_character_,
                     link = NA_character_, n = NA_integer_, k = NA_integer_, seed = NA_integer_) {
  rbind(df, data.frame(
    model     = as.character(model),
    trial     = as.integer(trial),
    converged = as.logical(converged),
    code      = as.integer(code),
    message   = as.character(message),
    link      = as.character(link),
    n         = as.integer(n),
    k         = as.integer(k),
    seed      = as.integer(seed),
    stringsAsFactors = FALSE
  ))
}

convergence_df <- data.frame(
  model=character(), trial=integer(), converged=logical(),
  code=integer(), message=character(),
  link=character(), n=integer(), k=integer(), seed=integer(),
  stringsAsFactors = FALSE
)



###############################################################
### MONTE CARLO LOOP
###############################################################
## (Put this once before the Monte Carlo loop)
convergence_df <- data.frame(
  model = character(), iter = integer(), ok = logical(),
  code = integer(), message = character(),
  link = character(), n = integer(), k = integer(), seed = integer(),
  stringsAsFactors = FALSE
)

## Helper to append one row (adjust if your add_conv differs)
add_conv <- function(df, model, iter, ok, code, message, link, n, k, seed) {
  r <- data.frame(model=model, iter=iter, ok=ok,
                  code=code, message=message, link=link,
                  n=n, k=k, seed=seed, stringsAsFactors=FALSE)
  rbind(df, r)
}



for (h in 1:trials) {
  
  GS <- sample_cal_beta(n, k, theta, link = link, cubic = FALSE)
  x  <- GS$x
  y  <- GS$Y
  y0 <- GS$Y0
  
  ######################################
  ### 1. LINEAR
  ######################################
  est <- tryCatch(Linear_cal(x, y, y0), error = function(e) e)
  if (inherits(est, "error")) {
    convergence_df <- add_conv(convergence_df, "Linear", h, FALSE, NA_integer_, conditionMessage(est), link, n, k, seedNum)
  } else {
    convergence_df <- add_conv(convergence_df, "Linear", h, TRUE, 0L, NA_character_, link, n, k, seedNum)
    linear.coeff[h]     <- est$x0
    linear.stderr[h]    <- est$SE.x0
    linear.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
    LL <- est$x0 - 1.96*est$SE.x0
    UL <- est$x0 + 1.96*est$SE.x0
    linear.coverage[h] <- (x0_true >= LL & x0_true <= UL)
  }
  
  ######################################
  ### 2. QUADRATIC
  ######################################
  est <- tryCatch(Quadratic_ML_cal(x, y, y0), error = function(e) e)
  if (inherits(est, "error")) {
    convergence_df <- add_conv(convergence_df, "Quadratic", h, FALSE, NA_integer_, conditionMessage(est), link, n, k, seedNum)
  } else {
    convergence_df <- add_conv(convergence_df, "Quadratic", h, TRUE, 0L, NA_character_, link, n, k, seedNum)
    quadratic.coeff[h]     <- est$x0
    quadratic.stderr[h]    <- est$SE.x0
    quadratic.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
    LL <- est$x0 - 1.96*est$SE.x0
    UL <- est$x0 + 1.96*est$SE.x0
    quadratic.coverage[h] <- (x0_true >= LL & x0_true <= UL)
  }
  
  ######################################
  ### 3. BCM (Beta_cal)
  ######################################
  est <- tryCatch(Beta_cal(x, y, y0, link = link),
                  error = function(e) e)
  
  if (inherits(est, "error")) {
    # Error path: record and move on
    convergence_df <- add_conv(
      convergence_df, "BCM", h, FALSE, NA_integer_,
      conditionMessage(est), link, n, k, seedNum
    )
    bcm.converge[h] <- NA_integer_
    
  } else {
    # Success path
    conv_code <- if (is.null(est$convergence)) NA_integer_ else as.integer(est$convergence)
    ok        <- isTRUE(conv_code == 0)
    msg       <- if (!is.null(est$message)) as.character(est$message) else NA_character_
    
    convergence_df <- add_conv(
      convergence_df, "BCM", h, ok, conv_code, msg, link, n, k, seedNum
    )
    bcm.converge[h] <- conv_code
    
    if (ok) {
      x0_hat <- est$estimates[3]
      se_x0  <- est$SE[3]
      
      bcm.coeff[h]     <- x0_hat
      bcm.stderr[h]    <- se_x0
      bcm.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
      
      LL <- x0_hat - 1.96 * se_x0
      UL <- x0_hat + 1.96 * se_x0
      bcm.coverage[h] <- (x0_true >= LL & x0_true <= UL)
    } else {
      # Converged != 0 (e.g., maxit): treat like non-estimated
      bcm.coverage[h] <- NA
    }
  }
  
  
  ######################################
  ### 4. WLS
  ######################################
  est <- tryCatch(WLS_cal(x, y, y0), error = function(e) e)
  if (inherits(est, "error")) {
    convergence_df <- add_conv(convergence_df, "WLS", h, FALSE, NA_integer_, conditionMessage(est), link, n, k, seedNum)
  } else {
    convergence_df <- add_conv(convergence_df, "WLS", h, TRUE, 0L, NA_character_, link, n, k, seedNum)
    wls.coeff[h]     <- est$x0
    wls.stderr[h]    <- est$SE.x0
    wls.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
    LL <- est$x0 - 1.96*est$SE.x0
    UL <- est$x0 + 1.96*est$SE.x0
    wls.coverage[h] <- (x0_true >= LL & x0_true <= UL)
  }
  
  ######################################
  ### 5. One-stage Beta Regression
  ######################################
  est <- tryCatch(Beta1stage_cal(x, y, y0, link = link), error = function(e) e)
  if (inherits(est, "error")) {
    convergence_df <- add_conv(convergence_df, "Beta1stage", h, FALSE, NA_integer_, conditionMessage(est), link, n, k, seedNum)
  } else {
    convergence_df <- add_conv(convergence_df, "Beta1stage", h, TRUE, 0L, NA_character_, link, n, k, seedNum)
    beta1.coeff[h]     <- est$x0
    beta1.stderr[h]    <- est$SE.x0
    beta1.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
    LL <- est$x0 - 1.96*est$SE.x0
    UL <- est$x0 + 1.96*est$SE.x0
    beta1.coverage[h] <- (x0_true >= LL & x0_true <= UL)
  }
  
  ######################################
  ### 6. Beta GAM
  ######################################
  est <- tryCatch(BetaGAM_cal(x, y, y0, B = 300), error = function(e) e)
  if (inherits(est, "error")) {
    convergence_df <- add_conv(convergence_df, "BetaGAM", h, FALSE, NA_integer_, conditionMessage(est), link, n, k, seedNum)
  } else {
    convergence_df <- add_conv(convergence_df, "BetaGAM", h, TRUE, 0L, NA_character_, link, n, k, seedNum)
    gam.coeff[h]     <- est$x0
    gam.stderr[h]    <- est$SE.x0
    gam.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
    LL <- est$x0 - 1.96*est$SE.x0
    UL <- est$x0 + 1.96*est$SE.x0
    gam.coverage[h] <- (x0_true >= LL & x0_true <= UL)
  }
  
  ######################################
  ### 7. Neural Network
  ######################################
  est <- tryCatch(Neural_cal(x, y, y0, B = 300), error = function(e) e)
  if (inherits(est, "error")) {
    convergence_df <- add_conv(convergence_df, "NeuralNet", h, FALSE, NA_integer_, conditionMessage(est), link, n, k, seedNum)
  } else {
    convergence_df <- add_conv(convergence_df, "NeuralNet", h, TRUE, 0L, NA_character_, link, n, k, seedNum)
    ann.coeff[h]     <- est$x0
    ann.stderr[h]    <- est$SE.x0
    ann.Criteria[h,] <- c(est$AIC, est$BIC, est$HQ)
    LL <- est$x0 - 1.96*est$SE.x0
    UL <- est$x0 + 1.96*est$SE.x0
    ann.coverage[h] <- (x0_true >= LL & x0_true <= UL)
  }
 print(h) 
} # end trials loop
# end Monte Carlo loop

###############################################################
### SAVE RESULTS
###############################################################

filename.est <- file.path(out_dir, paste0("All_estim_est_out", seedNum, ".csv"))
write.csv(
  cbind(
    linear.coeff, quadratic.coeff, bcm.coeff,
    wls.coeff, beta1.coeff, gam.coeff, ann.coeff
  ),
  file = filename.est,
  row.names = FALSE
)

filename.stderr <- file.path(out_dir, paste0("All_stderr_out", seedNum, ".csv"))
write.csv(
  cbind(
    linear.stderr, quadratic.stderr, bcm.stderr,
    wls.stderr, beta1.stderr, gam.stderr, ann.stderr
  ),
  file = filename.stderr,
  row.names = FALSE
)

filename.criteria <- file.path(out_dir, paste0("All_criteria_out", seedNum, ".csv"))
write.csv(
  cbind(
    linear.Criteria, quadratic.Criteria, bcm.Criteria,
    wls.Criteria, beta1.Criteria, gam.Criteria, ann.Criteria
  ),
  file = filename.criteria,
  row.names = FALSE
)

filename.coverage <- file.path(out_dir, paste0("All_coverage", seedNum, ".csv"))
write.csv(
  cbind(
    linear.coverage, quadratic.coverage, bcm.coverage,
    wls.coverage, beta1.coverage, gam.coverage, ann.coverage
  ),
  file = filename.coverage,
  row.names = FALSE
)

## unified convergence log (optional but nice)
filename.convergence <- file.path(out_dir, paste0("convergence_log", seedNum, ".csv"))
write.csv(convergence_df, file = filename.convergence, row.names = FALSE)
