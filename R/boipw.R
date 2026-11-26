## =====================================================================
## Bounded odd inverse Pareto Weibull (BOIPW) distribution on (0, 1)
##
##   G(y) = 1 - [ Kbar(y) / (beta + (1 - beta) Kbar(y)) ]^alpha,
##   K(y) = exp{-theta (-log y)^gamma},  Kbar = 1 - K,
##   alpha, beta, theta, gamma > 0.
##
## Functions follow the usual R naming convention:
##   dBOIPW (pdf), pBOIPW (cdf), qBOIPW (quantile), rBOIPW (random numbers),
##   hBOIPW (hazard rate), sBOIPW (survival), mBOIPW (raw moments).
## =====================================================================

## internal helper: t = -log y, w = K(y) = exp(-theta t^gamma), u = 1 - w
.boipw_tw <- function(y, theta, gamma) {
  t <- -log(y)
  w <- exp(-theta * t^gamma)
  u <- -expm1(-theta * t^gamma)
  list(t = t, w = w, u = u)
}

.check_par <- function(alpha, beta, theta, gamma) {
  if (any(c(alpha, beta, theta, gamma) <= 0))
    stop("all parameters must be positive")
}

## probability density function
dBOIPW <- function(x, alpha, beta, theta, gamma, log = FALSE) {
  .check_par(alpha, beta, theta, gamma)
  out <- rep(-Inf, length(x))
  ok <- x > 0 & x < 1
  y <- x[ok]
  z <- .boipw_tw(y, theta, gamma)
  ## log g(y), written in a numerically stable form
  lg <- log(alpha) + log(beta) + log(theta) + log(gamma) +
    (gamma - 1) * log(z$t) - theta * z$t^gamma - 2 * log(z$u) - log(y) -
    (alpha + 1) * log1p(beta * z$w / z$u)
  out[ok] <- lg
  if (log) out else exp(out)
}

## survival function S(y) = 1 - G(y)
sBOIPW <- function(q, alpha, beta, theta, gamma) {
  .check_par(alpha, beta, theta, gamma)
  out <- ifelse(q <= 0, 1, 0)
  ok <- q > 0 & q < 1
  z <- .boipw_tw(q[ok], theta, gamma)
  out[ok] <- exp(-alpha * log1p(beta * z$w / z$u))
  out
}

## cumulative distribution function
pBOIPW <- function(q, alpha, beta, theta, gamma, lower.tail = TRUE) {
  s <- sBOIPW(q, alpha, beta, theta, gamma)
  if (lower.tail) 1 - s else s
}

## hazard rate function h(y) = g(y) / S(y)
hBOIPW <- function(x, alpha, beta, theta, gamma) {
  .check_par(alpha, beta, theta, gamma)
  z <- .boipw_tw(x, theta, gamma)
  D <- beta + (1 - beta) * z$u
  alpha * beta * theta * gamma * z$t^(gamma - 1) * z$w / (x * z$u * D)
}

## quantile function (closed form, equation (8) of the paper)
qBOIPW <- function(p, alpha, beta, theta, gamma) {
  .check_par(alpha, beta, theta, gamma)
  if (any(p < 0 | p > 1)) stop("p must lie in [0, 1]")
  q <- (1 - p)^(1 / alpha)
  w <- beta * q / (1 - (1 - beta) * q)
  exp(-(-log1p(-w) / theta)^(1 / gamma))
}

## random number generation by the inverse transform method
rBOIPW <- function(n, alpha, beta, theta, gamma) {
  qBOIPW(runif(n), alpha, beta, theta, gamma)
}

## r-th raw moment E(Y^r) = int_0^1 Q(p)^r dp
mBOIPW <- function(r, alpha, beta, theta, gamma) {
  integrate(function(p) qBOIPW(p, alpha, beta, theta, gamma)^r, 0, 1,
            subdivisions = 1000L, rel.tol = 1e-10)$value
}

## mean, median, variance, skewness and kurtosis (Table 1 of the paper)
summaryBOIPW <- function(alpha, beta, theta, gamma) {
  m <- sapply(1:4, mBOIPW, alpha = alpha, beta = beta, theta = theta, gamma = gamma)
  mu <- m[1]; v <- m[2] - mu^2
  sk <- (m[3] - 3 * mu * m[2] + 2 * mu^3) / v^1.5
  ku <- (m[4] - 4 * mu * m[3] + 6 * mu^2 * m[2] - 3 * mu^4) / v^2
  c(mean = mu, median = qBOIPW(0.5, alpha, beta, theta, gamma),
    variance = v, skewness = sk, kurtosis = ku)
}

## ---------------------------------------------------------------------
## Sub-models and competitors used in the paper
##   BOIPE       : gamma = 1                       (alpha, beta, theta)
##   BMOEW       : alpha = 1                       (beta, theta, gamma)
##   UEW         : beta = 1                        (alpha, theta, gamma)
##                 unit exponentiated Weibull, G(y) = 1 - [1 - K(y)]^alpha
##   BMOEE       : alpha = gamma = 1               (beta, theta)
##   Kumaraswamy : beta = gamma = 1                (alpha, theta)
##                 G(y) = 1 - (1 - y^theta)^alpha
##   Unit-Weibull: alpha = beta = 1                (theta, gamma)
##                 G(y) = exp{-theta (-log y)^gamma}
##   EKw         : exponentiated Kumaraswamy (a, b, c), not nested in BOIPW
##                 G(y) = [1 - (1 - y^a)^b]^c
##   UBXII       : unit Burr-XII (a, b), not nested in BOIPW
##                 G(y) = [1 + (-log y)^a]^(-b)
##   Beta        : Beta(a, b), not nested in BOIPW
## ---------------------------------------------------------------------
MODELS <- list(
  BOIPW       = c("alpha", "beta", "theta", "gamma"),
  BOIPE       = c("alpha", "beta", "theta"),
  BMOEW       = c("beta", "theta", "gamma"),
  UEW         = c("alpha", "theta", "gamma"),
  BMOEE       = c("beta", "theta"),
  Kumaraswamy = c("alpha", "theta"),
  UW          = c("theta", "gamma"),
  EKw         = c("a", "b", "c"),
  UBXII       = c("a", "b"),
  Beta        = c("a", "b")
)
NON_NESTED <- c("EKw", "UBXII", "Beta")

## map the free parameters of a sub-model to (alpha, beta, theta, gamma)
full_par <- function(model, p) {
  switch(model,
    BOIPW       = p,
    BOIPE       = c(p[1], p[2], p[3], 1),
    BMOEW       = c(1, p[1], p[2], p[3]),
    UEW         = c(p[1], 1, p[2], p[3]),
    BMOEE       = c(1, p[1], p[2], 1),
    Kumaraswamy = c(p[1], 1, p[2], 1),
    UW          = c(1, 1, p[1], p[2]),
    stop("unknown model"))
}

model_logpdf <- function(model, y, p) {
  if (model == "Beta") return(dbeta(y, p[1], p[2], log = TRUE))
  if (model == "EKw") {
    ## numerically stable form: 1 - y^a = -expm1(a log y), v = 1 - (1 - y^a)^b
    l1 <- log(-expm1(p[1] * log(y)))          # log(1 - y^a)
    lv <- log(-expm1(p[2] * l1))              # log(1 - (1 - y^a)^b)
    return(log(p[1]) + log(p[2]) + log(p[3]) + (p[1] - 1) * log(y) + (p[2] - 1) * l1 +
             (p[3] - 1) * lv)
  }
  if (model == "UBXII") {
    t <- -log(y)
    return(log(p[1] * p[2]) - log(y) + (p[1] - 1) * log(t) - (p[2] + 1) * log1p(t^p[1]))
  }
  q <- full_par(model, p)
  dBOIPW(y, q[1], q[2], q[3], q[4], log = TRUE)
}

model_cdf <- function(model, y, p) {
  if (model == "Beta") return(pbeta(y, p[1], p[2]))
  if (model == "EKw") return(exp(p[3] * log(-expm1(p[2] * log(-expm1(p[1] * log(y)))))))
  if (model == "UBXII") return((1 + (-log(y))^p[1])^(-p[2]))
  q <- full_par(model, p)
  pBOIPW(y, q[1], q[2], q[3], q[4])
}

model_pdf <- function(model, y, p) exp(model_logpdf(model, y, p))

## random numbers from any model
model_random <- function(model, n, p) {
  u <- runif(n)
  if (model == "Beta") return(rbeta(n, p[1], p[2]))
  if (model == "EKw") return((1 - (1 - u^(1 / p[3]))^(1 / p[2]))^(1 / p[1]))
  if (model == "UBXII") return(exp(-(u^(-1 / p[2]) - 1)^(1 / p[1])))
  q <- full_par(model, p)
  qBOIPW(u, q[1], q[2], q[3], q[4])
}
