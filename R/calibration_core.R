# R/calibration_core.R
# All core functions used in simulations and applications

# libraries needed by core
suppressPackageStartupMessages({
  library(betareg)
  library(numDeriv)
  # mgcv, nnet are loaded inside functions via requireNamespace
})


# ============================================
# Data generation (BCM), supports linear or cubic predictor
# ============================================
sample_cal_beta <- function(n, k, theta,
                            link = c("logit", "probit", "cloglog"),
                            cubic = FALSE) {
  link <- match.arg(link)
  x <- seq(0, 2.5, length.out = n)
  
  if (!cubic) {
    stopifnot(length(theta) == 5L)
    b0 <- theta[1]; b1 <- theta[2]; x0 <- theta[3]; phi <- theta[4]; phi0 <- theta[5]
    eta <- b0 + b1 * x
    eta0 <- b0 + b1 * x0
  } else {
    stopifnot(length(theta) == 7L)
    b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]
    x0 <- theta[5]; phi <- theta[6]; phi0 <- theta[7]
    eta  <- b0 + b1*x + b2*x^2 + b3*x^3
    eta0 <- b0 + b1*x0 + b2*x0^2 + b3*x0^3
  }
  
  ginv <- switch(link,
                 logit   = function(u) 1/(1+exp(-u)),
                 probit  = function(u) pnorm(u),
                 cloglog = function(u) 1 - exp(-exp(u)))
  mu  <- .clamp01(ginv(eta))
  mu0 <- .clamp01(ginv(eta0))
  
  Y  <- rbeta(n, mu * phi, (1 - mu) * phi)
  Y0 <- rbeta(k, mu0 * phi0, (1 - mu0) * phi0)
  list(x = x, Y = Y, Y0 = Y0)
}


# ============================================
# MODELS
# ============================================

# ---- Two-stage BCM (proposed) ----
Beta_cal <- function(x, y, y0,
                     link = c("logit", "probit", "cloglog"),
                     start.v = NULL,
                     clamp_eps = .Machine$double.eps) {
  link <- match.arg(link)
  n1 <- length(y); n2 <- length(y0); n_total <- n1 + n2
  y  <- .clamp01(y,  clamp_eps); y0 <- .clamp01(y0, clamp_eps)
  
  g     <- switch(link, logit=qlogis, probit=qnorm, cloglog=function(mu) log(-log(1-mu)))
  ginv  <- switch(link, logit=function(eta) 1/(1+exp(-eta)),
                  probit=pnorm,
                  cloglog=function(eta) 1 - exp(-exp(eta)))
  
  if (is.null(start.v)) {
    fit <- betareg::betareg(y ~ x, link = link)
    a0  <- unname(fit$coefficients$mean[1]); b0 <- unname(fit$coefficients$mean[2])
    phi_s1 <- as.numeric(fit$coefficients$precision)
    y0bar <- mean(y0); y0bar <- .clamp01(y0bar, clamp_eps)
    x0.start <- (g(y0bar) - a0)/b0
    start.v <- c(a0, b0, x0.start, max(phi_s1,1), max(phi_s1,1))
  }
  
  loglik <- function(th) {
    alpha <- th[1]; beta <- th[2]; x0 <- th[3]
    phi <- th[4];  phi0 <- th[5]
    if (phi <= 0 || phi0 <= 0) return(-Inf)
    mu  <- .clamp01(ginv(alpha + beta * x))
    mu0 <- .clamp01(ginv(alpha + beta * x0))
    ll1 <- sum(lgamma(phi) - lgamma(phi*mu) - lgamma(phi*(1-mu)) +
                 (phi*mu - 1)*log(y) + (phi*(1-mu)-1)*log(1-y))
    ll2 <- sum(lgamma(phi0) - lgamma(phi0*mu0) - lgamma(phi0*(1-mu0)) +
                 (phi0*mu0 - 1)*log(y0) + (phi0*(1-mu0)-1)*log(1-y0))
    ll1 + ll2
  }
  
  est <- optim(start.v, loglik, method="L-BFGS-B",
               lower=c(-Inf,-Inf,-Inf,1e-8,1e-8),
               control=list(fnscale=-1, maxit=1000),
               hessian=TRUE)
  
  np <- length(start.v)
  the <- est$par; names(the) <- c("alpha","beta","x0","phi","phi0")
  
  SE <- rep(NA_real_, np)
  if (is.matrix(est$hessian) && all(is.finite(est$hessian))) {
    invH <- tryCatch(solve(-est$hessian), error=function(e) NULL)
    if (!is.null(invH)) SE <- sqrt(pmax(diag(invH), 0))
  }
  
  logL <- est$value
  AIC <- -2*logL + 2*np
  BIC <- -2*logL + np*log(n_total)
  HQ  <- -2*logL + 2*np*log(log(n_total))
  
  list(model="BCM", link=link,
       x0 = the["x0"], SE.x0 = SE[3],
       AIC=AIC, BIC=BIC, HQ=HQ, AICc = AIC + (2*np*(np+1))/(n_total - np - 1),
       ok = (est$convergence==0))
}

# ---- Linear calibration (EURACHEM baseline) ----
Linear_cal <- function(x, y, y0) {
  n <- length(y); k <- length(y0); N <- n + k
  fit <- lm(y ~ x)
  w <- coef(fit); a <- w[1]; b <- w[2]
  y0b <- mean(y0)
  x0  <- (y0b - a)/b
  # pooled variance across stages (classical two-stage LS)
  RSS <- sum((y - a - b*x)^2) + sum((y0 - y0b)^2)
  sigma2 <- RSS/(N - 3)  # (a,b, and x0 implied) – keeps same flavor you used
  Sxx <- sum( (x - mean(x))^2 )
  SE.x0 <- sqrt( sigma2 / b^2 * ( 1/k + 1/n + ((mean(x)-x0)^2)/Sxx ) )
  
  logL <- -(N/2)*log(2*pi) - (N/2)*log(sigma2) - RSS/(2*sigma2)
  p <- 4; # (a,b,x0,sigma2)
  AIC <- -2*logL + 2*p
  BIC <- -2*logL + p*log(N)
  HQ  <- -2*logL + 2*p*log(log(N))
  AICc <- AIC + (2*p*(p+1))/(N - p - 1)
  list(model="Linear", link="gaussian", x0=x0, SE.x0=SE.x0,
       AIC=AIC,BIC=BIC,HQ=HQ,AICc=AICc, ok=TRUE)
}

# ---- Quadratic calibration (EURACHEM alternative) ----
Quadratic_ML_cal <- function(x, y, y0) {
  n <- length(y); k <- length(y0); N <- n + k
  z <- x^2
  fit <- lm(y ~ x + z)
  w <- coef(fit)
  y0b <- mean(y0)
  
  # root for beta0 + beta1*x0 + beta2*x0^2 = y0_bar
  root_fun <- function(xx) w[1] + w[2]*xx + w[3]*xx^2 - y0b
  x0 <- tryCatch(uniroot(root_fun, c(0, max(10, max(x)*2)), tol=1e-12)$root,
                 error=function(e) NA_real_)
  if (is.na(x0)) {
    # fallback by shape
    disc <- (w[2]/(2*w[3]))^2 - (w[1]-y0b)/w[3]
    x0 <- if (isTRUE(all(diff(y) >= 0))) { # increasing
      -(w[2]/(2*w[3])) + sqrt(disc)
    } else {
      -(w[2]/(2*w[3])) - sqrt(disc)
    }
  }
  
  yhat <- w[1] + w[2]*x + w[3]*z
  sigma2_ml <- ( sum((y - yhat)^2) + sum((y0 - y0b)^2) ) / N
  
  # loglik at (w, x0, sigma2_ml)
  loglik <- function(theta) {
    b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; x0h <- theta[4]; s2 <- theta[5]
    if (s2 <= 0) return(-Inf)
    res1 <- sum((y - (b0 + b1*x + b2*x^2))^2)
    res2 <- sum((y0 - (b0 + b1*x0h + b2*x0h^2))^2)
    RSE  <- res1 + res2
    -(N/2)*log(2*pi) - (N/2)*log(s2) - RSE/(2*s2)
  }
  
  theta_hat <- c(w, x0, sigma2_ml)
  Hinv <- tryCatch(-solve(hessian(loglik, x=theta_hat)), error=function(e) NULL)
  SE.x0 <- if (is.null(Hinv)) NA_real_ else sqrt(max(Hinv[4,4], 0))
  
  logL <- loglik(theta_hat)
  p <- length(theta_hat)
  AIC <- -2*logL + 2*p
  BIC <- -2*logL + p*log(N)
  HQ  <- -2*logL + 2*p*log(log(N))
  AICc <- AIC + (2*p*(p+1))/(N - p - 1)
  
  list(model="Quadratic", link="gaussian", x0=x0, SE.x0=SE.x0,
       AIC=AIC,BIC=BIC,HQ=HQ,AICc=AICc, ok=TRUE)
}

# ---- WLS (traditional) ----
WLS_cal <- function(x, y, y0) {
  if (length(unique(y)) < 3) return(list(model="WLS", link="gaussian", x0=NA, SE.x0=NA, AIC=NA,BIC=NA,HQ=NA,AICc=NA, ok=FALSE))
  n <- length(y); k <- length(y0); N <- n + k
  ols <- lm(y ~ x)
  r   <- residuals(ols)
  w   <- 1/(r^2 + 1e-8)
  fit <- lm(y ~ x, weights=w)
  a <- coef(fit)[1]; b <- coef(fit)[2]
  y0b <- mean(y0)
  x0 <- (y0b - a)/b
  
  RSSw <- sum(w*(y - fitted(fit))^2)
  sigma2 <- RSSw/(n - 2)
  xw <- sum(w*x)/sum(w)
  Sxx <- sum(w*(x - xw)^2)
  SE.x0 <- sqrt( sigma2 / b^2 * (1/k + 1/n + ((xw - x0)^2)/Sxx) )
  
  logL <- -(n/2)*log(2*pi*sigma2) - 0.5*sum(w*(y - fitted(fit))^2/sigma2)
  p <- 2; AIC <- -2*logL + 2*p
  BIC <- -2*logL + p*log(N)
  HQ  <- -2*logL + 2*p*log(log(N))
  AICc <- AIC + (2*p*(p+1))/(N - p - 1)
  
  list(model="WLS", link="gaussian", x0=x0, SE.x0=SE.x0,
       AIC=AIC,BIC=BIC,HQ=HQ,AICc=AICc, ok=TRUE)
}

# ---- 1-stage Beta regression (Stage 1 only) ----
Beta1stage_cal <- function(x, y, y0, link = c("logit","probit","cloglog")) {
  link <- match.arg(link)
  n <- length(y); k <- length(y0); N <- n + k
  
  fit <- betareg::betareg(y ~ x, link=link)
  a <- fit$coefficients$mean[1]; b <- fit$coefficients$mean[2]
  Sigma <- vcov(fit)
  
  mu0 <- .clamp01(mean(y0))
  g <- switch(link, logit=qlogis, probit=qnorm, cloglog=function(m) log(-log(1-m)))
  x0 <- (g(mu0) - a)/b
  
  grad <- c(-1/b, -(g(mu0)-a)/(b^2), 0)
  SE.x0 <- sqrt(as.numeric(t(grad) %*% Sigma %*% grad))
  
  logL <- as.numeric(logLik(fit)); p <- 3
  AIC <- -2*logL + 2*p
  BIC <- -2*logL + p*log(N)
  HQ  <- -2*logL + 2*p*log(log(N))
  AICc <- AIC + (2*p*(p+1))/(N - p - 1)
  
  list(model=sprintf("Beta1(%s)", link), link=link,
       x0=x0, SE.x0=SE.x0,
       AIC=AIC,BIC=BIC,HQ=HQ,AICc=AICc, ok=TRUE)
}

# ---- Beta GAM (logit link), parametric bootstrap SE ----
BetaGAM_cal <- function(x, y, y0, B=500) {
  if (!.opt_load("mgcv")) {
    return(list(model="BetaGAM", link="logit", x0=NA, SE.x0=NA,
                AIC=NA,BIC=NA,HQ=NA,AICc=NA, ok=FALSE,
                note="mgcv not installed"))
  }
  n <- length(y); k <- length(y0); N <- n + k
  dat <- data.frame(x=x, y=y)
  gm <- mgcv::gam(y ~ s(x, bs="tp"), data=dat, family=mgcv::betar())
  
  y0b <- mean(y0)
  x_min <- min(x); x_max <- max(x)
  inv_fun <- function(xx) as.numeric(predict(gm, newdata=data.frame(x=xx), type="response")) - y0b
  x0 <- tryCatch(uniroot(inv_fun, c(x_min, x_max), tol=1e-6)$root, error=function(e) NA_real_)
  
  # Boot SE
  mu_hat <- as.numeric(predict(gm, type="response"))
  phi_hat <- gm$family$getTheta(TRUE)
  x0_boot <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    yb <- rbeta(n, mu_hat*phi_hat, (1-mu_hat)*phi_hat)
    gm_b <- tryCatch(mgcv::gam(yb ~ s(x, bs="tp"), family=mgcv::betar()), error=function(e) NULL)
    if (is.null(gm_b)) next
    inv_b <- function(xx) as.numeric(predict(gm_b, newdata=data.frame(x=xx), type="response")) - y0b
    x0_boot[b] <- tryCatch(uniroot(inv_b, c(x_min, x_max), tol=1e-6)$root, error=function(e) NA_real_)
  }
  SE.x0 <- if (sum(!is.na(x0_boot))>1) sd(x0_boot, na.rm=TRUE) else NA_real_
  
  # ICs (beta likelihood)
  logL <- as.numeric(logLik(gm))
  p <- sum(gm$edf) + 1
  AIC <- -2*logL + 2*p
  BIC <- -2*logL + p*log(N)
  HQ  <- -2*logL + 2*p*log(log(N))
  AICc <- AIC + (2*p*(p+1))/(N - p - 1)
  
  list(model="BetaGAM", link="logit", x0=x0, SE.x0=SE.x0,
       AIC=AIC,BIC=BIC,HQ=HQ,AICc=AICc, ok=TRUE)
}

# ---- Simple NN comparator (optional; needs nnet) ----
Neural_cal <- function(x, y, y0, B=500) {
  if (!.opt_load("nnet")) {
    return(list(model="ANN", link="sigmoid", x0=NA, SE.x0=NA,
                AIC=NA,BIC=NA,HQ=NA,AICc=NA, ok=FALSE,
                note="nnet not installed"))
  }
  n <- length(y); k <- length(y0); N <- n + k
  x_min <- min(x); x_max <- max(x); xr <- x_max - x_min
  xs <- if (xr>0) (x - x_min)/xr else rep(0.5,n)
  dat <- data.frame(xs=xs, y=y)
  nn <- nnet::nnet(y ~ xs, data=dat, size=8, linout=FALSE, trace=FALSE, maxit=1000)
  mu <- as.numeric(predict(nn, type="raw"))
  # crude MOM precision
  eps <- .Machine$double.eps
  muC <- .clamp01(mu, eps)
  yC  <- .clamp01(y, eps)
  var_y <- var(yC)
  mean_mu_var <- mean(muC*(1-muC))
  phi <- max((if (var_y>0) mean_mu_var/var_y - 1 else 100), 1)
  
  y0b <- mean(y0)
  inv <- function(xx) {
    xs0 <- if (xr>0) (xx - x_min)/xr else 0.5
    as.numeric(predict(nn, newdata=data.frame(xs=xs0), type="raw")) - y0b
  }
  x0 <- tryCatch(uniroot(inv, c(x_min,x_max), tol=1e-6)$root, error=function(e) NA_real_)
  
  # bootstrap SE
  x0b <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    yb <- rbeta(n, muC*phi, (1-muC)*phi)
    datb <- data.frame(xs=xs, y=yb)
    nnb <- tryCatch(nnet::nnet(y ~ xs, data=datb, size=8, linout=FALSE, trace=FALSE, maxit=1000),
                    error=function(e) NULL)
    if (is.null(nnb)) next
    invb <- function(xx) {
      xs0 <- if (xr>0) (xx - x_min)/xr else 0.5
      as.numeric(predict(nnb, newdata=data.frame(xs=xs0), type="raw")) - y0b
    }
    x0b[b] <- tryCatch(uniroot(invb, c(x_min,x_max), tol=1e-6)$root, error=function(e) NA_real_)
  }
  SE.x0 <- if (sum(!is.na(x0b))>1) sd(x0b, na.rm=TRUE) else NA_real_
  
  # Approx Beta log-likelihood for ICs
  logL <- sum(dbeta(yC, muC*phi, (1-muC)*phi, log=TRUE))
  p <- 26  # approx #params as discussed
  AIC <- -2*logL + 2*p
  BIC <- -2*logL + p*log(N)
  HQ  <- -2*logL + 2*p*log(log(N))
  AICc <- AIC + (2*p*(p+1))/(N - p - 1)
  
  list(model="ANN", link="sigmoid", x0=x0, SE.x0=SE.x0,
       AIC=AIC,BIC=BIC,HQ=HQ,AICc=AICc, ok=TRUE)
}


# ============================================
# ORCHESTRATOR
# ============================================
fit_all_models <- function(x, y, y0,
                           beta_links = c("logit"),
                           include_gam = TRUE,
                           include_ann = TRUE,
                           include_beta1 = TRUE,
                           include_bcm = TRUE,
                           include_linear = TRUE,
                           include_quadratic = TRUE,
                           include_wls = TRUE) {
  res <- list()
  
  # Linear / Quadratic / WLS
  if (include_linear)    res$Linear    <- tryCatch(Linear_cal(x,y,y0), error=function(e) list(model="Linear", ok=FALSE, err=e$message))
  if (include_quadratic) res$Quadratic <- tryCatch(Quadratic_ML_cal(x,y,y0), error=function(e) list(model="Quadratic", ok=FALSE, err=e$message))
  if (include_wls)       res$WLS       <- tryCatch(WLS_cal(x,y,y0), error=function(e) list(model="WLS", ok=FALSE, err=e$message))
  
  # Beta1 (Stage 1 only)
  if (include_beta1) {
    for (L in beta_links) {
      key <- paste0("Beta1_", L)
      res[[key]] <- tryCatch(Beta1stage_cal(x,y,y0, link=L),
                             error=function(e) list(model=key, ok=FALSE, err=e$message))
    }
  }
  
  # BCM (two-stage)
  if (include_bcm) {
    for (L in beta_links) {
      key <- paste0("BCM_", L)
      res[[key]] <- tryCatch(Beta_cal(x,y,y0, link=L),
                             error=function(e) list(model=key, ok=FALSE, err=e$message))
    }
  }
  
  # GAM and ANN
  if (include_gam) res$BetaGAM <- tryCatch(BetaGAM_cal(x,y,y0, B=300),
                                           error=function(e) list(model="BetaGAM", ok=FALSE, err=e$message))
  if (include_ann) res$ANN     <- tryCatch(Neural_cal(x,y,y0, B=300),
                                           error=function(e) list(model="ANN", ok=FALSE, err=e$message))
  
  # Tidy table
  tidy <- do.call(rbind, lapply(res, function(z) {
    nm <- if (!is.null(z$model)) z$model else deparse(substitute(z))
    cbind.data.frame(
      Model = nm,
      x0    = if (!is.null(z$x0)) z$x0 else NA_real_,
      SE    = if (!is.null(z$SE.x0)) z$SE.x0 else NA_real_,
      LCL   = if (!is.null(z$x0) && !is.null(z$SE.x0)) z$x0 - .z975*z$SE.x0 else NA_real_,
      UCL   = if (!is.null(z$x0) && !is.null(z$SE.x0)) z$x0 + .z975*z$SE.x0 else NA_real_,
      AIC   = if (!is.null(z$AIC)) z$AIC else NA_real_,
      AICc  = if (!is.null(z$AICc)) z$AICc else NA_real_,
      BIC   = if (!is.null(z$BIC)) z$BIC else NA_real_,
      HQ    = if (!is.null(z$HQ)) z$HQ else NA_real_,
      ok    = if (!is.null(z$ok)) z$ok else FALSE,
      note  = if (!is.null(z$note)) z$note else NA_character_,
      row.names = NULL
    )
  }))
  
  # identify IC winners (among ok fits)
  ic_cols <- c("AIC","AICc","BIC","HQ")
  for (ic in ic_cols) {
    m <- tidy$ok & is.finite(tidy[[ic]])
    tidy[[paste0("best_",ic)]] <- FALSE
    if (any(m)) {
      idx <- which.min(tidy[[ic]][m])
      w   <- which(m)[idx]
      tidy[[paste0("best_",ic)]][w] <- TRUE
    }
  }
  tidy[order(tidy$Model), ]
}


