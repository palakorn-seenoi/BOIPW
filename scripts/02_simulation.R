## =====================================================================
## Monte Carlo simulation study (Table 2 of the paper)
## Run from the repository root:  Rscript scripts/02_simulation.R [R] [cores]
##   R     : number of replications (default 1000, as in the paper)
##   cores : number of CPU cores for parallel::mclapply (default 1;
##           use 1 on Windows, where forking is not available)
## The results depend on the random numbers, so they agree with Table 2
## up to Monte Carlo error but not to the last digit (the paper's run
## used a different random number generator).
## =====================================================================
source("R/boipw.R"); source("R/fit.R")
dir.create("results", showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000L
cores <- if (length(args) >= 2) as.integer(args[2]) else 1L

## log-likelihood and analytic score (equations (11)-(15)), on the log scale
nll_grad <- function(lp, y) {
  p <- exp(lp); a <- p[1]; b <- p[2]; th <- p[3]; ga <- p[4]; n <- length(y)
  t <- -log(y); tg <- t^ga; w <- exp(-th * tg); u <- -expm1(-th * tg)
  D <- 1 - (1 - b) * w; lt <- log(t)
  ll <- n * log(a * b * th * ga) + (ga - 1) * sum(lt) - th * sum(tg) +
    (a - 1) * sum(log(u)) + sum(t) - (a + 1) * sum(log(D))
  s <- c(n / a + sum(log(u)) - sum(log(D)),
         n / b - (a + 1) * sum(w / D),
         n / th - sum(tg) + (a - 1) * sum(tg * w / u) - (a + 1) * (1 - b) * sum(tg * w / D),
         n / ga + sum(lt) - th * sum(tg * lt) + (a - 1) * th * sum(tg * lt * w / u) -
           (a + 1) * (1 - b) * th * sum(tg * lt * w / D))
  if (!is.finite(ll) || any(!is.finite(s))) return(list(f = 1e10, g = rep(0, 4)))
  list(f = -ll, g = -s * p)
}

## The log-likelihood can have several local maxima, so each sample is
## fitted from 13 starting points: the true values, two random perturbations
## of them, all parameters equal to 1, and nine random points that do not
## use the true values. The estimate is the best of the 13 solutions.
## The function also returns the best log-likelihood found from the three
## true-value starts only and from the first eight starts only, and the best
## solution from the ten starts that do not use the true values (sensitivity
## analysis of Tables S1-S2 of the Supporting Information).
fit_sim <- function(y, true) {
  starts <- c(list(log(true)), replicate(2, log(true) + rnorm(4, 0, 0.5), simplify = FALSE),
              list(rep(0, 4)), replicate(9, rnorm(4, 0, 1), simplify = FALSE))
  val <- rep(Inf, length(starts)); par <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    r <- try(optim(starts[[i]], function(l) nll_grad(l, y)$f, function(l) nll_grad(l, y)$g,
                   method = "L-BFGS-B", lower = -10, upper = 10,
                   control = list(maxit = 2000)), silent = TRUE)
    if (!inherits(r, "try-error")) { val[i] <- r$value; par[[i]] <- r$par }
  }
  b <- which.min(val); b0 <- 3 + which.min(val[4:13])
  list(est = exp(par[[b]]), ll = -val[b], ll.true = -min(val[1:3]), ll.first8 = -min(val[1:8]),
       est.nt = exp(par[[b0]]))
}

settings <- list(S1 = c(0.9, 3.8, 0.8, 1.0), S2 = c(1.2, 1.2, 2.0, 1.5), S3 = c(1.6, 1.6, 1.4, 0.8))
ns <- c(30, 50, 100, 200, 300, 500)
pars <- c("alpha", "beta", "theta", "gamma")
out <- NULL; extra <- NULL
for (si in names(settings)) {
  true <- settings[[si]]
  for (n in ns) {
    set.seed(10000 * match(si, names(settings)) + n)
    seeds <- sample.int(.Machine$integer.max, R)
    one <- function(seed) {
      set.seed(seed)
      y <- rBOIPW(n, true[1], true[2], true[3], true[4])
      y <- pmin(pmax(y, 1e-12), 1 - 1e-12)
      f <- fit_sim(y, true); est <- f$est
      ## boundary hit: an estimate at the edge of the search box [e^-10, e^10]
      hit <- any(abs(log(est)) > 10 - 1e-3)
      ## 95% Wald intervals from the analytic observed information (Appendix B)
      se <- tryCatch(sqrt(diag(solve(-boipw_hessian(est, y)))),
                     error = function(e) rep(NA, 4), warning = function(w) rep(NA, 4))
      cover <- if (hit || any(!is.finite(se))) rep(NA, 4) else abs(est - true) <= 1.959964 * se
      ## Wald intervals on the log scale: log(est) +/- 1.96 se / est
      cover.log <- if (hit || any(!is.finite(se))) rep(NA, 4) else abs(log(est) - log(true)) <= 1.959964 * se / est
      ## distance between the fitted and the true cdf, sup_y |G(y; est) - G(y; true)|
      yg <- seq(0.0005, 0.9995, length.out = 2000)
      dist <- max(abs(pBOIPW(yg, est[1], est[2], est[3], est[4]) - pBOIPW(yg, true[1], true[2], true[3], true[4])))
      c(est, hit, cover, f$ll - f$ll.true, f$ll - f$ll.first8, dist, cover.log, f$est.nt)
    }
    A <- if (cores > 1) do.call(rbind, parallel::mclapply(seeds, one, mc.cores = cores))
         else t(sapply(seeds, one))
    E <- A[, 1:4]
    res <- data.frame(setting = si, n = n, parameter = pars, true = true,
                      bias = colMeans(E) - true,
                      median.bias = apply(E, 2, median) - true,
                      rmse = sqrt(colMeans(sweep(E, 2, true)^2)))
    E0 <- A[, 17:20]
    res$median.bias.nt <- apply(E0, 2, median) - true
    res$rmse.nt <- sqrt(colMeans(sweep(E0, 2, true)^2))
    print(res, digits = 4, row.names = FALSE)
    out <- rbind(out, res)
    ## Tables S1-S2 (Supporting Information): boundary hits, coverage of 95% Wald intervals (among the
    ## replications with an interior estimate and a positive definite
    ## information matrix), and the percentage of samples in which the
    ## true-value starts alone, or the first eight starts alone, would have
    ## given a log-likelihood lower by more than 0.01
    cp <- colMeans(A[, 6:9, drop = FALSE], na.rm = TRUE)
    cpl <- colMeans(A[, 13:16, drop = FALSE], na.rm = TRUE)
    extra <- rbind(extra, data.frame(setting = si, n = n,
      boundary.pct = 100 * mean(A[, 5]), wald.valid.pct = 100 * mean(!is.na(A[, 6])),
      cp.alpha = cp[1], cp.beta = cp[2], cp.theta = cp[3], cp.gamma = cp[4],
      true.start.lower.pct = 100 * mean(A[, 10] > 0.01),
      first8.lower.pct = 100 * mean(A[, 11] > 0.01),
      mean.cdf.distance = mean(A[, 12]),
      cplog.alpha = cpl[1], cplog.beta = cpl[2], cplog.theta = cpl[3], cplog.gamma = cpl[4],
      nt.same.pct = 100 * mean(apply(abs(log(A[, 17:20]) - log(A[, 1:4])), 1, max) < 1e-3)))
    print(tail(extra, 1), digits = 3, row.names = FALSE)
  }
}
write.csv(out, "results/table2_simulation.csv", row.names = FALSE)
write.csv(extra, "results/tableS1_S2_simulation.csv", row.names = FALSE)
cat("Table 2 written to results/table2_simulation.csv\n")
