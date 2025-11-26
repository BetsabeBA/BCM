# Beta Calibration Model (BCM)

This repository contains the full R implementation, simulation code, and application scripts for:

**Blas Achic & Tavares Cavalcante (2025)**  
*Beta Calibration for Bounded Signals: A Two-Stage Likelihood Framework for Analytical Chemistry*  
(Chemometrics and Intelligent Laboratory Systems — submitted)

The goal of this repository is **complete reproducibility** of all numerical results, simulation experiments, and real-data analyses presented in the manuscript.

---

## 📌 Overview

The Beta Calibration Model (BCM) is a **two-stage likelihood framework** for calibration problems where the analytical signal is a **bounded, heteroscedastic** measurement such as:

- spectrophotometric transmittance  
- ratio measurements  
- proportions or absorbance scaled to (0,1)

This repository includes:

### **Core BCM implementation**
- Two-stage beta likelihood  
- Logit, probit, and complementary log–log links  
- Fisher-based standard errors using the observed Hessian

### **Comparator Models**
| Category | Methods |
|---------|---------|
| Classical calibration | Linear, Quadratic, Weighted Least Squares (WLS) |
| Likelihood-based | One-stage Beta regression |
| Semi-parametric | Beta GAM (penalized splines) |
| Machine learning | Neural network (single-hidden-layer ANN) |

### **Simulation Framework**
- Monte Carlo simulation engine for **logit / probit / cloglog** links  
- Supports **linear** and **cubic** data-generating mechanisms  
- Implements **7 models** (BCM + 6 competitors)  
- Stores convergence logs, AIC/BIC/HQ model-selection statistics, and performance metrics (Bias, SE, STD, RMSE, CP)

### **Real-data Applications**
- Iron in water (transmittance)
- Tryptophan spectrophotometry (absorbance → transmittance transformation)

---

## 📁 Folder Structure

```text
BCM/
├─ README.md                     # This file
├─ LICENSE                       # MIT license
│
├─ R/
│  ├─ calibration_core.R         # All model implementations:
│  │                               BCM, Linear, Quadratic, WLS,
│  │                               Beta1-stage, BetaGAM, ANN
│  └─ helpers/                   # Optional future expansion
│
├─ simulation/
│  ├─ BCM_simulation.R           # Main Monte Carlo driver (500+ replicates)
│  ├─ BCM_summary.R              # Parses All_* outputs and produces tables
│  ├─ exec.sh                    # Linux batch runner (SLURM/cluster-ready)
│  ├─ log/                       # Runtime logs (contains .gitkeep)
│  ├─ n=5-logit/                 # Results (auto-generated)
│  ├─ n=5-probit/
│  ├─ n=5-cloglog/
│  ├─ n=10-*/                    # Etc. (all auto-generated)
│  └─ .gitkeep                   # Keeps simulation directory structure
│
└─ examples/
   ├─ BCM_applications.R         # Iron + Tryptophan case studies
   └─ data/                      # Optional raw data folder
      └─ .gitkeep
