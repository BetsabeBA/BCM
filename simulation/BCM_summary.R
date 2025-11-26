############################################################
## BCM_summary.R
## Summary for BCM simulation study
## - Reads All_*.csv outputs (7 models)
## - Uses per-model convergence info (if available)
## - Masks non-converged fits
## - Reports Bias, SE, SD, RMSE, CP and % best by AIC/BIC/HQ
############################################################

# rm(list = ls())

suppressPackageStartupMessages({
  library(xtable)
})

## ------------------ USER SETTINGS -------------------------

## Sample sizes n to summarize
n_grid   <-5      #c(5, 10, 200, 500)   # adjust as needed, or use a single value

## Link used in the simulation folder names: n=5-logit, n=5-probit, n=5-cloglog
link     <- "probit"             # "logit", "probit", or "cloglog"

## Base directory where all n=...-link folders live (relative to repo root)
base_dir <- "simulation"

## Seeds: 1 for k=3, 2 for k=10 (matching BCM_simulation.R)
seed_k3  <- 1
seed_k10 <- 2

## True x0 used in the simulation
x0_true  <- 1.25

## Model names in output columns (7 models)
model_names <- c("Linear","Quadratic","BCM","WLS","Beta1","BetaGAM","ANN")
n_models    <- length(model_names)

## ------------------ UTILITIES -----------------------------

## Coerce any vector to logical "OK?" flag (TRUE = converged)
is_ok_vec <- function(v) {
  if (is.logical(v)) return(ifelse(is.na(v), FALSE, v))
  vn <- suppressWarnings(as.numeric(v))
  if (any(!is.na(vn))) return(ifelse(is.na(vn), FALSE, vn %in% c(0,1)))
  txt <- tolower(trimws(as.character(v)))
  okset <- c("0","ok","success","true","converged")
  ifelse(is.na(txt), FALSE, txt %in% okset)
}

## Read CSV and drop an index column if present (X, Row, etc.)
read_drop_index <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  df <- read.csv(path, header = TRUE, check.names = FALSE)
  if (ncol(df) >= 2 && names(df)[1] %in% c("X","X.1","row","row.names","Row","V1","index",".id",".row")) {
    df <- df[, -1, drop = FALSE]
  }
  df
}

## Ensure data.frame has exactly 'need' columns (pad with NA or truncate)
ensure_ncols <- function(df, need, prefix = "V") {
  have <- ncol(df)
  if (have == need) return(df)
  if (have > need)  return(df[, seq_len(need), drop = FALSE])
  add <- need - have
  pad <- as.data.frame(matrix(NA, nrow = nrow(df), ncol = add))
  names(pad) <- paste0(prefix, seq_len(add))
  cbind(df, pad)
}

## Read unified per-model convergence log (if available)
## Expected columns: model, iter/rep/h, ok/convergence/status
## Returns: matrix [n_trials x n_models] logical (TRUE = OK)
read_unified_convergence <- function(path, n_trials, model_names) {
  if (!file.exists(path)) return(NULL)
  df <- tryCatch(read.csv(path, header = TRUE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  
  pick <- function(cands) {
    ix <- which(tolower(names(df)) %in% cands)
    if (length(ix)) ix[1] else NA_integer_
  }
  
  col_model <- pick(c("model"))
  col_rep   <- pick(c("rep","h","iter","trial","i"))
  col_ok    <- pick(c("ok","convergence","status","conv_ok","conv"))
  if (is.na(col_model) || is.na(col_rep) || is.na(col_ok)) return(NULL)
  
  mdl <- as.character(df[[col_model]])
  rep <- suppressWarnings(as.integer(df[[col_rep]]))
  okv <- is_ok_vec(df[[col_ok]])
  
  M <- matrix(TRUE, nrow = n_trials, ncol = length(model_names))
  colnames(M) <- model_names
  
  rep[is.na(rep) | rep < 1 | rep > n_trials] <- NA
  
  for (j in seq_along(mdl)) {
    if (is.na(rep[j])) next
    m <- mdl[j]
    if (!(m %in% model_names)) next
    M[rep[j], m] <- okv[j]
  }
  M
}

## Legacy BCM-only convergence file (if unified one not present)
read_bcm_convergence_vec <- function(path, n_trials) {
  if (!file.exists(path)) return(rep(TRUE, n_trials))
  df <- tryCatch(read.csv(path, header = TRUE, check.names = FALSE), error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(rep(TRUE, n_trials))
  pick_cols <- c("bcm.converge","convergence","converged","status","ok","success",
                 "BCM","bcm","Beta","beta","conv")
  col_idx <- which(tolower(names(df)) %in% tolower(pick_cols))
  v <- if (length(col_idx)) df[[col_idx[1]]] else df[[1]]
  ok <- is_ok_vec(v)
  if (length(ok) != n_trials) ok <- rep(ok, length.out = n_trials)
  ok
}

## Read all result files for a given (nsize, seed, link)
read_results_for_seed <- function(nsize, seed, base_dir, model_names, link) {
  dir_n <- file.path(base_dir, sprintf("n=%d-%s", nsize, link))
  
  est_fp  <- file.path(dir_n, paste0("All_estim_est_out", seed, ".csv"))
  se_fp   <- file.path(dir_n, paste0("All_stderr_out",   seed, ".csv"))
  crit_fp <- file.path(dir_n, paste0("All_criteria_out", seed, ".csv"))
  cov_fp  <- file.path(dir_n, paste0("All_coverage",     seed, ".csv"))
  
  x0_est <- read_drop_index(est_fp)
  se_est <- read_drop_index(se_fp)
  crit   <- read_drop_index(crit_fp)
  cover  <- read_drop_index(cov_fp)
  
  ## 7 models => 7 estimate/SE/coverage cols; 21 criteria cols
  x0_est <- ensure_ncols(x0_est, n_models, prefix = "x0_pad")
  se_est <- ensure_ncols(se_est, n_models, prefix = "se_pad")
  cover  <- ensure_ncols(cover,  n_models, prefix = "cov_pad")
  crit   <- ensure_ncols(crit,   3 * n_models, prefix = "crit_pad")
  
  colnames(x0_est) <- model_names
  colnames(se_est) <- model_names
  colnames(cover)  <- model_names
  
  crit_names <- as.vector(unlist(lapply(model_names, function(m) paste0(m, c(".aic",".bic",".hq")))))
  if (ncol(crit) == length(crit_names)) {
    names(crit) <- crit_names
  } else {
    nm <- names(crit)
    for (m in model_names) {
      for (suf in c(".aic",".bic",".hq")) {
        tgt <- paste0(m, suf)
        hit <- which(tolower(nm) == tolower(tgt))
        if (length(hit) == 0) {
          hit <- which(grepl(paste0("^", m, ".*", sub("^\\.", "\\\\.", suf), "$"),
                             nm, ignore.case = TRUE))
        }
        if (length(hit) == 1) nm[hit] <- tgt
      }
    }
    names(crit) <- nm
  }
  
  n_trials <- nrow(x0_est)
  
  ## Convergence logs
  conv_log_fp <- file.path(dir_n, paste0("convergence_log", seed, ".csv"))  # unified
  conv_old_fp <- file.path(dir_n, paste0("convergence",     seed, ".csv"))  # legacy
  
  conv_ok_mat <- read_unified_convergence(conv_log_fp, n_trials, model_names)
  if (is.null(conv_ok_mat)) {
    conv_ok_mat <- matrix(TRUE, nrow = n_trials, ncol = n_models)
    colnames(conv_ok_mat) <- model_names
    bcm_vec <- read_bcm_convergence_vec(conv_old_fp, n_trials)
    conv_ok_mat[, "BCM"] <- bcm_vec
  }
  
  list(
    x0_est   = x0_est,
    se_est   = se_est,
    crit     = crit,
    cover    = cover,
    conv_ok  = conv_ok_mat,
    n_trials = n_trials,
    dir      = dir_n,
    seed     = seed
  )
}

## Summary metrics for x0 by model
summarize_x0 <- function(x0_est, se_est, cover, x0_true, model_names) {
  bias    <- apply(x0_est - x0_true, 2, mean, na.rm = TRUE)
  sd_x0   <- apply(x0_est - x0_true, 2, sd,   na.rm = TRUE)
  mean_se <- apply(se_est,           2, mean, na.rm = TRUE)
  rmse    <- sqrt(apply((x0_est - x0_true)^2, 2, mean, na.rm = TRUE))
  cp      <- apply(cover, 2, mean, na.rm = TRUE) * 100
  
  out <- rbind(bias, mean_se, sd_x0, rmse, cp)
  rownames(out) <- c("Bias", "SE", "STD", "RMSE", "CP")
  colnames(out) <- model_names
  out
}

## Helper: argmin with NA handling
which_min_ignore_na <- function(v) {
  w <- replace(v, is.na(v), Inf)
  if (all(is.infinite(w))) return(NA_integer_)
  which.min(w)
}

## Model-selection % (AIC/BIC/HQ) across ALL models
selection_perc <- function(crit_mat, model_names) {
  get_mat <- function(suffix) {
    cols <- paste0(model_names, ".", suffix)
    out  <- matrix(NA_real_, nrow = nrow(crit_mat), ncol = length(cols))
    colnames(out) <- model_names
    for (j in seq_along(cols)) {
      if (cols[j] %in% names(crit_mat)) out[, j] <- crit_mat[[cols[j]]]
    }
    out
  }
  AIC_mat <- get_mat("aic")
  BIC_mat <- get_mat("bic")
  HQ_mat  <- get_mat("hq")
  
  pick <- function(M) {
    ix <- apply(M, 1, which_min_ignore_na)
    nm <- ifelse(is.na(ix), NA_character_, model_names[ix])
    sapply(model_names, function(m) mean(nm == m, na.rm = TRUE) * 100)
  }
  
  out <- rbind("AIC %" = pick(AIC_mat),
               "BIC %" = pick(BIC_mat),
               "HQ  %" = pick(HQ_mat))
  colnames(out) <- model_names
  out
}

## Apply a logical mask (FALSE => NA) to one model column
mask_col <- function(x, model, mask) {
  x[[model]][!mask] <- NA
  x
}

## Mask non-converged rows for EVERY model, across x0_est, se_est, cover, crit
mask_nonconvergence <- function(res, model_names) {
  for (m in model_names) {
    ok <- res$conv_ok[, m]
    res$x0_est <- mask_col(res$x0_est, m, ok)
    res$se_est <- mask_col(res$se_est, m, ok)
    res$cover  <- mask_col(res$cover,  m, ok)
    for (suf in c(".aic",".bic",".hq")) {
      cn <- paste0(m, suf)
      if (cn %in% names(res$crit)) res$crit[[cn]][!ok] <- NA
    }
  }
  res
}

## Count non-convergence per model
count_nonconv <- function(conv_ok_mat, model_names) {
  colSums(!conv_ok_mat)
}

## ------------------ MAIN LOOP -----------------------------

out_all <- list()

for (nsize in n_grid) {
  
  ## ---- k = 3 ----
  res_k3 <- read_results_for_seed(nsize, seed_k3, base_dir, model_names, link)
  non_k3 <- count_nonconv(res_k3$conv_ok, model_names)
  res_k3 <- mask_nonconvergence(res_k3, model_names)
  summary_k3 <- summarize_x0(res_k3$x0_est, res_k3$se_est, res_k3$cover, x0_true, model_names)
  select_k3  <- selection_perc(res_k3$crit, model_names)
  
  ## ---- k = 10 ----
  res_k10 <- read_results_for_seed(nsize, seed_k10, base_dir, model_names, link)
  non_k10 <- count_nonconv(res_k10$conv_ok, model_names)
  res_k10 <- mask_nonconvergence(res_k10, model_names)
  summary_k10 <- summarize_x0(res_k10$x0_est, res_k10$se_est, res_k10$cover, x0_true, model_names)
  select_k10  <- selection_perc(res_k10$crit, model_names)
  
  out_all[[paste0("n=", nsize)]] <- list(
    summary = list(k3 = summary_k3,  k10 = summary_k10),
    select  = list(k3 = select_k3,    k10 = select_k10),
    nonconv = list(
      k3  = data.frame(Model = model_names,
                       NonConverged = as.integer(non_k3),
                       Total = res_k3$n_trials,
                       Pct = 100 * as.integer(non_k3) / res_k3$n_trials),
      k10 = data.frame(Model = model_names,
                       NonConverged = as.integer(non_k10),
                       Total = res_k10$n_trials,
                       Pct = 100 * as.integer(non_k10) / res_k10$n_trials)
    )
  )
}

## ------------------ PRINT RESULTS -------------------------

for (nm in names(out_all)) {
  cat("\n==================== ", nm, " (link = ", link, ") ====================\n", sep = "")
  
  cat("\n-- Non-convergence per model (k = 3) --\n")
  print(out_all[[nm]]$nonconv$k3, row.names = FALSE)
  cat("\n-- Non-convergence per model (k = 10) --\n")
  print(out_all[[nm]]$nonconv$k10, row.names = FALSE)
  
  cat("\n-- % best by criteria (k = 3) --\n")
  print(round(out_all[[nm]]$select$k3, 1))
  cat("\n-- % best by criteria (k = 10) --\n")
  print(round(out_all[[nm]]$select$k10, 1))
  
  cat("\n-- Summary (Bias, SE, STD, RMSE, CP) (k = 3) --\n")
  print(round(out_all[[nm]]$summary$k3, 3))
  cat("\n-- Summary (Bias, SE, STD, RMSE, CP) (k = 10) --\n")
  print(round(out_all[[nm]]$summary$k10, 3))
  
  cat("\n-- LaTeX (Summary k=3) --\n")
  print(xtable(out_all[[nm]]$summary$k3, digits = 3), include.rownames = TRUE)
  cat("\n-- LaTeX (Summary k=10) --\n")
  print(xtable(out_all[[nm]]$summary$k10, digits = 3), include.rownames = TRUE)
}

## ------------------ Optional save of one block -------------------------
## ex <- out_all[["n=5"]]
## write.csv(ex$nonconv$k3,  file.path(base_dir, "nonconv_k3_n5.csv"),  row.names = FALSE)
## write.csv(ex$nonconv$k10, file.path(base_dir, "nonconv_k10_n5.csv"), row.names = FALSE)
