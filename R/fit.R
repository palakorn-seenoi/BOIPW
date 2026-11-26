## =====================================================================
## Maximum likelihood estimation, goodness-of-fit and likelihood ratio
## tests for the BOIPW distribution and its competitors.
## Requires R/boipw.R. Base R only.
## =====================================================================

## negative log-likelihood on the log-parameter scale (keeps parameters > 0)
.nll_log <- function(lp, model, y) {
  v <- -sum(model_logpdf(model, y, exp(lp)))
  if (is.finite(v)) v else 1e10
}

## Fit a model by maximum likelihood with a multi-start strategy:
## each start is refined by Nelder-Mead followed by L-BFGS-B; all
## log-parameters are kept in [-12, 12] (parameters in [e^-12, e^12]),
## and the best solution is kept. A parameter estimate at the edge of this
## box indicates that the likelihood is maximised on the boundary of the
## parameter space (a limiting case of the model).
fit_model <- function(model, y, nstart = 60, seed = 1) {
  k <- length(MODELS[[model]])
  set.seed(seed)
  starts <- c(list(rep(0, k)), replicate(nstart, rnorm(k, 0, 1.5), simplify = FALSE))
  best <- NULL
  for (s in starts) {
    r1 <- try(optim(s, .nll_log, model = model, y = y, method = "Nelder-Mead",
                    control = list(maxit = 20000, reltol = 1e-14)), silent = TRUE)
    if (inherits(r1, "try-error")) next
    r2 <- try(optim(pmin(pmax(r1$par, -12), 12), .nll_log, model = model, y = y,
                    method = "L-BFGS-B", lower = -12, upper = 12,
                    control = list(maxit = 5000, factr = 1)), silent = TRUE)
    ## the Nelder-Mead result is accepted only if it stays inside the search box
    ok1 <- all(abs(r1$par) <= 12)
    r <- if (!inherits(r2, "try-error") && (!ok1 || r2$value < r1$value)) r2 else if (ok1) r1 else next
    if (is.null(best) || r$value < best$value - 1e-9) best <- r
  }
  par <- exp(best$par); names(par) <- MODELS[[model]]
  at_bound <- any(abs(best$par) > 12 - 1e-3)
  nll <- best$value; n <- length(y)
  ## standard errors from the observed information J = -Hessian, computed
  ## from the analytic second derivatives (Appendix B of the paper)
  ## (NA when the MLE lies on the boundary of the parameter space, where the
  ## information matrix is singular or cannot be evaluated)
  ## the information matrix is inverted on the log-parameter scale,
  ## J_log = diag(par) J diag(par), which is algebraically equivalent but much
  ## better conditioned when the parameters differ by several orders of
  ## magnitude; se(par) = par * se(log par)
  se <- tryCatch({ Jl <- -diag(par, k) %*% loglik_hessian(model, par, y) %*% diag(par, k)
                   v <- diag(solve(Jl)); if (any(v <= 0)) stop("not positive definite"); par * sqrt(v) },
                 error = function(e) rep(NA, k), warning = function(w) rep(NA, k))
  names(se) <- names(par)
  if (at_bound) se[] <- NA
  list(model = model, par = par, se = se, loglik = -nll, at_bound = at_bound,
       AIC = 2 * nll + 2 * k, BIC = 2 * nll + k * log(n), n = n)
}

## ---------- Anderson-Darling p-value (Marsaglia and Marsaglia, 2004) ----------
.adinf <- function(z) {
  if (z < 2) {
    exp(-1.2337141 / z) / sqrt(z) *
      (2.00012 + (0.247105 - (0.0649821 - (0.0347962 - (0.011672 - 0.00168691 * z) * z) * z) * z) * z)
  } else {
    exp(-exp(1.0776 - (2.30695 - (0.43424 - (0.082433 - (0.008056 - 0.0003146 * z) * z) * z) * z) * z))
  }
}
.errfix <- function(n, x) {
  c <- 0.01265 + 0.1757 / n
  if (x < c) {
    t <- x / c; t <- sqrt(t) * (1 - t) * (49 * t - 102)
    return(t * (0.0037 / n^3 + 0.00078 / n^2 + 0.00006 / n))
  } else if (x < 0.8) {
    t <- (x - c) / (0.8 - c)
    t <- -0.00022633 + (6.54034 - (14.6538 - (14.458 - (8.259 - 1.91864 * t) * t) * t) * t) * t
    return(t * (0.04213 / n + 0.01365 / n^2))
  }
  (-130.2137 + (745.2337 - (1705.091 - (1950.646 - (1116.360 - 255.7844 * x) * x) * x) * x) * x) / n
}
ad_test <- function(u) {
  u <- sort(pmin(pmax(u, 1e-300), 1 - 1e-16)); n <- length(u); i <- seq_len(n)
  A <- -n - mean((2 * i - 1) * (log(u) + log1p(-rev(u))))
  x <- .adinf(A)
  c(statistic = A, p.value = 1 - (x + .errfix(n, x)))
}

## ---------- Cramer-von Mises p-value (Csorgo and Faraway, 1996) ----------
.cvm_inf <- function(x) {
  tot <- 0; k <- 0
  repeat {
    u <- exp(lgamma(k + 0.5) - lgamma(k + 1)) / (pi^1.5 * sqrt(x))
    y <- 4 * k + 1; q <- y^2 / (16 * x)
    z <- u * sqrt(y) * exp(-q) * besselK(q, 0.25)
    tot <- tot + z; k <- k + 1
    if (abs(z) < 1e-7) break
  }
  tot
}
.psi1_mod <- function(x) {
  ed2 <- function(y) { z <- y^2 / 4; exp(-z) * (y / 2)^1.5 * (besselK(z, 1/4) + besselK(z, 3/4)) / sqrt(pi) }
  ed3 <- function(y) { z <- y^2 / 4; exp(-z) / sqrt(pi) * (y / 2)^2.5 *
                         (2 * besselK(z, 1/4) + 3 * besselK(z, 3/4) - besselK(z, 5/4)) }
  tot <- 0; k <- 0
  repeat {
    m <- 2 * k + 1; sx <- 2 * sqrt(x); y1 <- x^0.75; y2 <- x^1.25
    g1 <- gamma(k + 0.5); g3 <- gamma(k + 1.5)
    A <- m * g1 * ed2((4 * k + 3) / sx) / (9 * y1) +
      g1 * ed3((4 * k + 1) / sx) / (72 * y2) +
      2 * (m + 2) * g3 * ed3((4 * k + 5) / sx) / (12 * y2) +
      7 * m * g1 * ed2((4 * k + 1) / sx) / (144 * y1) +
      7 * m * g1 * ed2((4 * k + 5) / sx) / (144 * y1)
    z <- -A / (pi * gamma(k + 1))
    tot <- tot + z; k <- k + 1
    if (abs(z) < 1e-7) break
  }
  tot
}
cvm_test <- function(u) {
  u <- sort(u); n <- length(u); i <- seq_len(n)
  W <- 1 / (12 * n) + sum(((2 * i - 1) / (2 * n) - u)^2)
  cdf <- if (W <= 1 / (12 * n)) 0 else if (W >= n / 3) 1 else
    .cvm_inf(W) * (1 + 1 / (12 * n)) + .psi1_mod(W) / n
  c(statistic = W, p.value = max(0, 1 - cdf))
}

## K-S, A-D and Cramer-von Mises statistics with p-values
gof_stats <- function(fit, y) {
  u <- model_cdf(fit$model, y, fit$par)
  ## exact two-sided p-value (as in the paper); data I contains a few tied values,
  ## for which ks.test() would otherwise switch to the asymptotic approximation
  ks <- suppressWarnings(ks.test(u, "punif", exact = TRUE))
  ad <- ad_test(u); cv <- cvm_test(u)
  c(KS = unname(ks$statistic), KS.p = ks$p.value,
    AD = unname(ad["statistic"]), AD.p = unname(ad["p.value"]),
    W2 = unname(cv["statistic"]), W2.p = unname(cv["p.value"]))
}

## likelihood ratio test of a nested sub-model against the BOIPW fit
LR_DF <- c(BOIPE = 1, BMOEW = 1, UEW = 1, BMOEE = 2, Kumaraswamy = 2, UW = 2)
LR_H0 <- c(BOIPE = "gamma = 1", BMOEW = "alpha = 1", UEW = "beta = 1", BMOEE = "alpha = gamma = 1",
           Kumaraswamy = "beta = gamma = 1", UW = "alpha = beta = 1")
lr_test <- function(fit_full, fit_sub) {
  L <- max(0, 2 * (fit_full$loglik - fit_sub$loglik))
  df <- LR_DF[[fit_sub$model]]
  c(Lambda = L, df = df, p.value = pchisq(L, df, lower.tail = FALSE))
}

## ---------- analytic score and Hessian of the BOIPW log-likelihood ----------
## (equations (20)-(23) and Appendix B, (B1)-(B8), of the paper)
boipw_score <- function(par, y) {
  a <- par[1]; b <- par[2]; th <- par[3]; ga <- par[4]; n <- length(y)
  t <- -log(y); L <- log(t); s <- t^ga; K <- exp(-th * s); Kb <- -expm1(-th * s)
  D <- 1 - (1 - b) * K; A <- K / Kb; C <- K / D
  c(alpha = n / a + sum(log(Kb)) - sum(log(D)),
    beta  = n / b - (a + 1) * sum(C),
    theta = n / th - sum(s) + (a - 1) * sum(s * A) - (a + 1) * (1 - b) * sum(s * C),
    gamma = n / ga + sum(L) - th * sum(s * L) + (a - 1) * th * sum(s * L * A) -
            (a + 1) * (1 - b) * th * sum(s * L * C))
}
boipw_hessian <- function(par, y) {
  a <- par[1]; b <- par[2]; th <- par[3]; ga <- par[4]; n <- length(y)
  t <- -log(y); L <- log(t); s <- t^ga; K <- exp(-th * s); Kb <- -expm1(-th * s)
  D <- 1 - (1 - b) * K; A <- K / Kb; C <- K / D
  H <- matrix(0, 4, 4, dimnames = list(names(par), names(par)))
  H[1, 1] <- -n / a^2
  H[1, 2] <- -sum(C)
  H[1, 3] <- sum(s * (A - (1 - b) * C))
  H[1, 4] <- th * sum(s * L * (A - (1 - b) * C))
  H[2, 2] <- -n / b^2 + (a + 1) * sum(C^2)
  H[2, 3] <- (a + 1) * sum(s * C / D)
  H[2, 4] <- (a + 1) * th * sum(s * L * C / D)
  H[3, 3] <- -n / th^2 - (a - 1) * sum(s^2 * A / Kb) + (a + 1) * (1 - b) * sum(s^2 * C / D)
  H[3, 4] <- -sum(s * L) + (a - 1) * sum(s * L * A * (1 - th * s / Kb)) -
             (a + 1) * (1 - b) * sum(s * L * C * (1 - th * s / D))
  H[4, 4] <- -n / ga^2 - th * sum(s * L^2) + (a - 1) * th * sum(s * L^2 * A * (1 - th * s / Kb)) -
             (a + 1) * (1 - b) * th * sum(s * L^2 * C * (1 - th * s / D))
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  H
}

## Hessian of the log-likelihood of any model used in the paper (analytic
## except for the non-nested EKw and unit Burr-XII competitors):
## nested sub-models use the corresponding sub-matrix of the BOIPW Hessian
## (the fixed parameters are held at 1); Beta(a, b) uses trigamma functions.
loglik_hessian <- function(model, par, y) {
  n <- length(y)
  if (model == "Beta") {
    a <- par[1]; b <- par[2]; tab <- trigamma(a + b)
    return(n * matrix(c(tab - trigamma(a), tab, tab, tab - trigamma(b)), 2, 2))
  }
  if (model %in% c("EKw", "UBXII"))      # no closed form used: numerical Hessian
    ## computed on the log scale (keeps the parameters positive) and
    ## transformed back; exact at the MLE, where the score is zero
    return(-optimHess(log(par), function(l) -sum(model_logpdf(model, y, exp(l))),
                      control = list(ndeps = rep(1e-4, length(par)))) / outer(par, par))
  H <- boipw_hessian(full_par(model, par), y)
  idx <- match(MODELS[[model]], c("alpha", "beta", "theta", "gamma"))
  H[idx, idx, drop = FALSE]
}
