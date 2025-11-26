## =======================================================
## BCM_applications.R
## Applications of BCM and comparators to real datasets:
## - Iron in water (transmittance)
## - Tryptophan (Samples A and B, transmittance)
## Requires: R/calibration_apply.R in the same repo
## =======================================================

## Adjust this to where you cloned the repo
## setwd("path/to/BCM")

source("R/calibration_apply.R")

## -------------------------------------------------------
## 1) Iron in water – transmittance (as in the paper)
## -------------------------------------------------------

iron_data <- matrix(
  c(0.2, 0.1351, 0.73266,
    1.0, 0.7169, 0.19191,
    1.5, 1.0846, 0.08230,
    2.0, 1.4416, 0.03617,
    2.5, 1.6849, 0.02066),
  nrow = 5, ncol = 3, byrow = TRUE
)
colnames(iron_data) <- c("x", "absorbance", "transmittance")

x_iron  <- iron_data[, "x"]
y_iron  <- iron_data[, "transmittance"]

x0_true_iron <- 0.2
y0_abs_trans <- c(0.1519, 0.70486)  # (absorbance, transmittance)
y0_iron      <- y0_abs_trans[2]     # transmittance

cat("\n==================================================\n")
cat("IRON DATASET – TRANSMITTANCE\n")
cat("==================================================\n")

for (L in c("logit", "probit", "cloglog")) {
  cat("\n--- Link:", L, "---\n")
  res_iron <- fit_all_models(
    x  = x_iron,
    y  = y_iron,
    y0 = y0_iron,
    beta_links   = L,
    include_gam  = TRUE,
    include_ann  = TRUE
  )
  print(res_iron, digits = 4)
}

## -------------------------------------------------------
## 2) Tryptophan – Sample A (absorvance → transmittance)
##    true x0 = 7
## -------------------------------------------------------

tryA_data <- matrix(
  c(0.0, 0.0356,
    5.0, 0.3068,
    # 7.0, 0.3440,
    10.0, 0.3980,
    14.0, 0.3670,
    15.0, 0.3860,
    20.0, 0.6020,
    25.0, 0.6680,
    29.0, 0.8470),
  nrow = 8, ncol = 2, byrow = TRUE
)
colnames(tryA_data) <- c("x", "absorbance")

x_A  <- tryA_data[, "x"]
y_A  <- tryA_data[, "absorbance"]
y0_A <- 0.3440   # absorbance, true x0 = 7

# transform to transmittance
y_A_t  <- 10^(2 - y_A)/100
y0_A_t <- 10^(2 - y0_A)/100

cat("\n==================================================\n")
cat("TRYPTOPHAN – SAMPLE A (TRANSMITTANCE)\n")
cat("==================================================\n")

for (L in c("logit", "probit", "cloglog")) {
  cat("\n--- Link:", L, "---\n")
  res_A <- fit_all_models(
    x  = x_A,
    y  = y_A_t,
    y0 = y0_A_t,
    beta_links   = L,
    include_gam  = TRUE,
    include_ann  = TRUE
  )
  print(res_A, digits = 4)
}

## -------------------------------------------------------
## 3) Tryptophan – Sample B (absorvance → transmittance)
##    true x0 = 14
## -------------------------------------------------------

tryB_data <- matrix(
  c(0.0, 0.0356,
    5.0, 0.3068,
    7.0, 0.3440,
    10.0, 0.3980,
    #14.0, 0.3670,
    15.0, 0.3860,
    20.0, 0.6020,
    25.0, 0.6680,
    29.0, 0.8470),
  nrow = 8, ncol = 2, byrow = TRUE
)
colnames(tryB_data) <- c("x", "absorbance")

x_B  <- tryB_data[, "x"]
y_B  <- tryB_data[, "absorbance"]
y0_B <- 0.3670   # absorbance, true x0 = 14

# transform to transmittance
y_B_t  <- 10^(2 - y_B)/100
y0_B_t <- 10^(2 - y0_B)/100

cat("\n==================================================\n")
cat("TRYPTOPHAN – SAMPLE B (TRANSMITTANCE)\n")
cat("==================================================\n")

for (L in c("logit", "probit", "cloglog")) {
  cat("\n--- Link:", L, "---\n")
  res_B <- fit_all_models(
    x  = x_B,
    y  = y_B_t,
    y0 = y0_B_t,
    beta_links   = L,
    include_gam  = TRUE,
    include_ann  = TRUE
  )
  print(res_B, digits = 4)
}

## End of file
