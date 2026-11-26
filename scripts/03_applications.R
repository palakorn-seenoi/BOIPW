## =====================================================================
## Applications to real data (Tables 3-10 and Figures 5-6 of the paper)
## Run from the repository root:  Rscript scripts/03_applications.R
## =====================================================================
source("R/boipw.R"); source("R/fit.R")
dir.create("results", showWarnings = FALSE)

models <- c("BOIPW", "BOIPE", "BMOEW", "UEW", "BMOEE", "Kumaraswamy", "UW", "EKw", "UBXII", "Beta")
datasets <- list(
  I  = list(name = "ncyoungmale", y = read.csv("data/ncyoungmale.csv")$proportion),
  II = list(name = "catholic", y = read.csv("data/catholic.csv")$proportion)
)

describe <- function(y) {
  n <- length(y); m <- mean(y); s <- sd(y)
  ## bias-corrected sample skewness and excess kurtosis (as in scipy.stats)
  g1 <- sqrt(n * (n - 1)) / (n - 2) * (mean((y - m)^3) / mean((y - m)^2)^1.5)
  g2 <- (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * (mean((y - m)^4) / mean((y - m)^2)^2 - 3) + 6)
  q <- quantile(y, c(0.25, 0.5, 0.75), type = 7)
  c(n = n, min = min(y), Q1 = q[[1]], median = q[[2]], mean = m, Q3 = q[[3]],
    max = max(y), sd = s, skewness = g1, excess.kurtosis = g2)
}

all_fits <- list()
for (lab in names(datasets)) {
  y <- datasets[[lab]]$y; nm <- datasets[[lab]]$name
  cat("\n================ Data", lab, "(", nm, ", n =", length(y), ") ================\n")

  ## descriptive statistics
  d <- describe(y); print(round(d, 4))
  write.csv(data.frame(statistic = names(d), value = d),
            sprintf("results/table_desc_data%s.csv", lab), row.names = FALSE)

  ## maximum likelihood fits
  fits <- lapply(models, fit_model, y = y, nstart = 60)
  names(fits) <- models
  all_fits[[lab]] <- fits

  ## MLEs and standard errors
  est <- do.call(rbind, lapply(fits, function(f) {
    out <- setNames(rep(NA, 8), c("alpha", "beta", "theta", "gamma",
                                  "se.alpha", "se.beta", "se.theta", "se.gamma"))
    nmz <- names(f$par)
    if (f$model %in% NON_NESTED)   # report (a, b, c) in the alpha, beta, theta columns
      nmz <- c("alpha", "beta", "theta")[seq_along(nmz)]
    out[nmz] <- f$par; out[paste0("se.", nmz)] <- f$se; out }))
  cat("\nMLEs (standard errors):\n"); print(signif(est, 5))
  write.csv(data.frame(model = models, est), sprintf("results/table_mle_data%s.csv", lab), row.names = FALSE)

  ## log-likelihood, AIC, BIC and goodness-of-fit statistics
  gof <- do.call(rbind, lapply(fits, function(f)
    c(loglik = f$loglik, AIC = f$AIC, BIC = f$BIC, gof_stats(f, y))))
  cat("\nGoodness of fit:\n"); print(round(gof, 4))
  write.csv(data.frame(model = models, gof), sprintf("results/table_gof_data%s.csv", lab), row.names = FALSE)

  ## likelihood ratio tests against BOIPW
  lr <- do.call(rbind, lapply(names(LR_DF), function(m) lr_test(fits$BOIPW, fits[[m]])))
  lr <- data.frame(submodel = names(LR_DF), H0 = LR_H0[names(LR_DF)], lr)
  cat("\nLikelihood ratio tests:\n"); print(lr, digits = 4, row.names = FALSE)
  write.csv(lr, sprintf("results/table_lr_data%s.csv", lab), row.names = FALSE)

  ## Figure: (a) histogram with fitted densities, (b) empirical cdf with fitted cdfs
  ## BMOEE, EKw and UBXII are omitted for clarity
  cols <- c(BOIPW = "#c00000", UEW = "#8c564b", BMOEW = "#1f4e79", BOIPE = "#2e7d32",
            Kumaraswamy = "#ef6c00", Beta = "#7b1fa2", UW = "#455a64")
  ltys <- c(BOIPW = 1, UEW = 3, BMOEW = 2, BOIPE = 4, Kumaraswamy = 3, Beta = 5, UW = 6)
  lo <- 0; hi <- if (lab == "I") 0.3 else 1
  xx <- seq(1e-4, hi - 1e-4, length.out = 1000)
  png(sprintf("results/figure_data%s.png", lab), width = 3.6, height = 4.8, units = "in", res = 300)
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), cex = 0.7)
  dens <- sapply(names(cols), function(m) model_pdf(m, xx, fits[[m]]$par))
  hist(y, breaks = if (lab == "I") seq(0.06, 0.26, by = 0.01) else 20, freq = FALSE, col = "grey85", border = "grey40",
       xlim = c(lo, hi), ylim = c(0, if (lab == "I") 60 else 25),
       main = sprintf("Data %s: histogram and fitted densities", lab), xlab = "y")
  for (m in names(cols)) lines(xx, dens[, m], col = cols[m], lty = ltys[m], lwd = 1.3)
  legend("topright", legend = names(cols), col = cols, lty = ltys, bty = "n", cex = 0.8)
  mtext("(a)", side = 3, adj = -0.12, line = 0.5, font = 2)
  plot(ecdf(y), do.points = FALSE, verticals = TRUE, xlim = c(lo, hi),
       main = sprintf("Data %s: empirical and fitted cdfs", lab), xlab = "y", ylab = "cumulative probability")
  for (m in names(cols)) lines(xx, model_cdf(m, xx, fits[[m]]$par), col = cols[m], lty = ltys[m], lwd = 1.3)
  mtext("(b)", side = 3, adj = -0.12, line = 0.5, font = 2)
  dev.off()
}
saveRDS(all_fits, "results/fits.rds")
cat("\nResults written to the results/ folder.\n")
