## =====================================================================
## Numerical verification of the analytical results of the paper
## Run from the repository root:  Rscript scripts/00_verify_formulas.R
##
## Each check compares a closed-form expression with an independent
## numerical computation (numerical integration, numerical
## differentiation or Monte Carlo) and prints the largest error over
## several parameter settings. Base R only.
## =====================================================================
source("R/boipw.R"); source("R/fit.R")

settings <- list(c(1.2, 1.2, 2, 1.5), c(0.5, 0.5, 0.5, 0.5), c(1.6, 1.6, 1.4, 0.8),
                 c(3, 0.2, 0.5, 2), c(2, 5, 1, 0.5), c(0.9, 1.7, 0.8, 1.0))
res <- data.frame(check = character(), max_error = numeric(), tolerance = numeric())
report <- function(check, err, tol) {
  res[nrow(res) + 1, ] <<- list(check, err, tol)
}
integ <- function(f, lo = 0, hi = 1) {
  r <- tryCatch(integrate(f, lo, hi, rel.tol = 1e-10, abs.tol = 0, subdivisions = 5000L)$value,
                error = function(e) NA)
  if (is.na(r)) r <- integrate(f, lo, hi, rel.tol = 1e-8, subdivisions = 5000L)$value
  r
}
## pdf with respect to u = K(y): g(y) dy = gu(u) du, gu(u) = a b (1-u)^(a-1) / (1-(1-b)u)^(a+1).
## Integrals over y are computed on the u scale, where the integrand is smoother.
QK <- function(u, th, ga) exp(-(-log(u) / th)^(1 / ga))
gu <- function(u, a, b) a * b * (1 - u)^(a - 1) / (1 - (1 - b) * u)^(a + 1)
yy <- c(0.03, 0.3, 0.7, 0.95)

## Proposition 1, eqs (4)-(7): pdf = dG/dy, h = g/S, total mass 1
e1 <- e2 <- e3 <- 0
for (s in settings) {
  G <- function(y) pBOIPW(y, s[1], s[2], s[3], s[4])
  g <- dBOIPW(yy, s[1], s[2], s[3], s[4])
  hstep <- 1e-6 * yy
  e1 <- max(e1, abs(g - (G(yy + hstep) - G(yy - hstep)) / (2 * hstep)) / g)
  e2 <- max(e2, abs(hBOIPW(yy, s[1], s[2], s[3], s[4]) - g / sBOIPW(yy, s[1], s[2], s[3], s[4])) / g)
  e3 <- max(e3, abs(integ(function(u) gu(u, s[1], s[2])) - 1))
}
report("Prop. 1: pdf (5) = derivative of cdf (4)  [rel.]", e1, 1e-6)
report("Prop. 1: hazard (7) = g/S                 [rel.]", e2, 1e-8)
report("Remark 1: integral of the pdf = 1", e3, 1e-8)

## Proposition 2: boundary behaviour (ratio to the asymptotic form -> 1)
e1 <- e2 <- 0
for (s in settings) {
  a <- s[1]; b <- s[2]; th <- s[3]; ga <- s[4]
  y <- 1 - 1e-9
  e1 <- max(e1, abs(dBOIPW(y, a, b, th, ga) / (a * ga * (th / b)^a * (1 - y)^(a * ga - 1)) - 1),
            abs(hBOIPW(y, a, b, th, ga) / (a * ga / (1 - y)) - 1))
  y <- 1e-300; t <- -log(y)
  ## compared on the log scale to avoid overflow: log g - log(asymptotic form) -> 0
  e2 <- max(e2, abs(dBOIPW(y, a, b, th, ga, log = TRUE) -
                      (log(a * b * th * ga) + (ga - 1) * log(t) + t - th * t^ga)))
}
report("Prop. 2(i): asymptotic form as y -> 1    [rel.]", e1, 1e-4)
report("Prop. 2(ii): asymptotic form as y -> 0   [log]", e2, 1e-3)

## Proposition 3, eq (8): G(Q(p)) = p
e1 <- 0
for (s in settings) {
  p <- c(0.001, 0.1, 0.5, 0.9, 0.999)
  e1 <- max(e1, abs(pBOIPW(qBOIPW(p, s[1], s[2], s[3], s[4]), s[1], s[2], s[3], s[4]) - p))
}
report("Prop. 3: quantile (8), G(Q(p)) = p", e1, 1e-8)

## Proposition 4, eq (9): three expressions for E(Y^r)
FX <- function(x, a, b, th, ga) { u <- -expm1(-th * x^ga); (u / (b + (1 - b) * u))^a }
e1 <- 0
for (s in settings) for (r in c(1, 2, 3.5)) {
  m1 <- integ(function(u) QK(u, s[3], s[4])^r * gu(u, s[1], s[2]))
  m2 <- integ(function(p) qBOIPW(p, s[1], s[2], s[3], s[4])^r)
  m3 <- r * integ(function(x) exp(-r * x) * FX(x, s[1], s[2], s[3], s[4]), 0, Inf)
  e1 <- max(e1, abs(m1 - m2), abs(m1 - m3))
}
report("Prop. 4: moment expressions (9)", e1, 1e-7)

## Proposition 5, eqs (10)-(13): mixture representations and moments
K  <- function(y, th, ga) exp(-th * (-log(y))^ga)
kf <- function(y, th, ga) th * ga * (-log(y))^(ga - 1) * K(y, th, ga) / y
J <- 0:400
e1 <- e2 <- e3 <- 0
for (s in settings) {
  a <- s[1]; b <- s[2]; th <- s[3]; ga <- s[4]
  g <- dBOIPW(yy, a, b, th, ga)
  if (b < 2) {                                  # (10)
    ser <- sapply(yy, function(y) sum(b * (1 - b)^J * kf(y, th, ga) * K(y, th, ga)^J *
                                        (1 - K(y, th, ga))^(a - 1) / beta(J + 1, a)))
    e1 <- max(e1, abs(ser - g) / g)
    m12 <- sum(sapply(0:150, function(j) b * (1 - b)^j / beta(j + 1, a) *
                        integ(function(u) QK(u, th, ga) * u^j * (1 - u)^(a - 1))))
    e3 <- max(e3, abs(m12 - integ(function(u) QK(u, th, ga) * gu(u, a, b))))
  }
  if (b > 0.5) {                                # (11)
    rho <- 1 - 1 / b
    v <- exp(lgamma(a + J) - lgamma(a) - lfactorial(J)) * b^(-a) * rho^J
    ser <- sapply(yy, function(y) sum(v * (a + J) * kf(y, th, ga) * (1 - K(y, th, ga))^(a + J - 1)))
    e2 <- max(e2, abs(ser - g) / g)
    m13 <- sum(sapply(0:150, function(j) v[j + 1] * (a + j) *
                        integ(function(u) QK(u, th, ga) * (1 - u)^(a + j - 1))))
    e3 <- max(e3, abs(m13 - integ(function(u) QK(u, th, ga) * gu(u, a, b))))
  }
}
report("Prop. 5: geometric expansion (10)        [rel.]", e1, 1e-8)
report("Prop. 5: negative binomial expansion (11) [rel.]", e2, 1e-8)
report("Moments from the expansions (12)-(13), r = 1", e3, 1e-6)

## Proposition 6, eq (16): R = alpha2 / (alpha1 + alpha2)
e1 <- e2 <- 0
set.seed(1)
for (s in settings) for (a2 in c(0.4, 2.3)) {
  a1 <- s[1]; b <- s[2]; th <- s[3]; ga <- s[4]
  ## integrate over u = K(y) to avoid the singularities of the pdf
  f <- function(u) gu(u, a1, b) * pBOIPW(QK(u, th, ga), a2, b, th, ga)
  e1 <- max(e1, abs(integ(f) - a2 / (a1 + a2)))
  N <- 2e5
  mc <- mean(rBOIPW(N, a2, b, th, ga) < rBOIPW(N, a1, b, th, ga))
  e2 <- max(e2, abs(mc - a2 / (a1 + a2)) / sqrt(0.25 / N))
}
report("Prop. 6: R by numerical integration", e1, 1e-7)
report("Prop. 6: R by simulation (in standard errors)", e2, 4)

## Score (20)-(23) and Hessian (B1)-(B8): central finite differences
dat <- list(read.csv("data/ncyoungmale.csv")$proportion, read.csv("data/catholic.csv")$proportion)
ll <- function(p, y) sum(dBOIPW(y, p[1], p[2], p[3], p[4], log = TRUE))
e1 <- e2 <- 0
for (y in dat) for (p in list(c(1.2, 1.2, 2, 1.5), c(0.9, 1.7, 0.8, 1.0), c(2, 0.6, 1.1, 1.3))) {
  names(p) <- c("alpha", "beta", "theta", "gamma")
  h <- 1e-5 * pmax(1, abs(p))
  num_s <- sapply(1:4, function(j) { e <- replace(numeric(4), j, h[j])
    (ll(p + e, y) - ll(p - e, y)) / (2 * h[j]) })
  num_H <- sapply(1:4, function(j) { e <- replace(numeric(4), j, h[j])
    (boipw_score(p + e, y) - boipw_score(p - e, y)) / (2 * h[j]) })
  sc <- boipw_score(p, y); H <- boipw_hessian(p, y)
  e1 <- max(e1, abs(sc - num_s) / pmax(1, abs(sc)))
  e2 <- max(e2, abs(H - num_H) / pmax(1, abs(H)))
}
report("Score functions (20)-(23)               [rel.]", e1, 1e-6)
report("Hessian (B1)-(B8)                       [rel.]", e2, 1e-6)

res$status <- ifelse(res$max_error <= res$tolerance, "OK", "CHECK")
res$max_error <- signif(res$max_error, 3)
print(res, row.names = FALSE, right = FALSE)
dir.create("results", showWarnings = FALSE)
write.csv(res, "results/verification.csv", row.names = FALSE)
if (any(res$status != "OK")) stop("some checks exceeded their tolerance") else
  cat("\nAll checks passed.\n")
