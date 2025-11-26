# Beta Calibration Model (BCM)

This repository contains the R code and example datasets used in the paper

> Blas Achic & Tavares Cavalcante, "Beta Calibration for Bounded Signals: A Two-Stage Likelihood Framework for Analytical Chemistry" (Chemometrics and Intelligent Laboratory Systems, submitted).

## Contents

- `R/calibration_apply.R`  
  Core implementation of the Beta Calibration Model (BCM) and all comparator
  methods:
  - Linear and quadratic calibration
  - Weighted least squares (WLS)
  - One-stage Beta regression calibration
  - Two-stage Beta Calibration Model (BCM)
  - Beta GAM comparator (`mgcv::gam` with `betar()`)
  - Neural network comparator (`nnet`)

  The main user-facing function is:

  ```r
  fit_all_models(x, y, y0, ...)
