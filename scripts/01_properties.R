## =====================================================================
## Distributional properties (Table 1 and Figures 2-4 of the paper)
## Run from the repository root:  Rscript scripts/01_properties.R
## =====================================================================
source("R/boipw.R")
dir.create("results", showWarnings = FALSE)

## ---- Table 1: mean, median, variance, skewness and kurtosis ----
settings <- rbind(c(0.5, 0.5, 0.5, 0.5), c(0.5, 8, 2, 1.2), c(1.2, 1.2, 2, 1.5),
                  c(1.6, 1.6, 1.4, 0.8), c(2, 5, 1, 0.5), c(3, 0.2, 0.5, 2),
                  c(0.9, 3.8, 0.8, 1.0))
tab1 <- t(apply(settings, 1, function(s) summaryBOIPW(s[1], s[2], s[3], s[4])))
tab1 <- data.frame(alpha = settings[, 1], beta = settings[, 2], theta = settings[, 3],
                   gamma = settings[, 4], round(tab1, 4))
print(tab1, row.names = FALSE)
write.csv(tab1, "results/table1_moments.csv", row.names = FALSE)

## check: the numerical moments agree with the pdf (mean by direct integration)
s <- settings[3, ]
cat("mean check:", integrate(function(y) y * dBOIPW(y, s[1], s[2], s[3], s[4]), 0, 1)$value,
    "vs", tab1$mean[3], "\n")

lab <- function(s) bquote(alpha == .(s[1]) * "," ~ beta == .(s[2]) * "," ~
                            theta == .(s[3]) * "," ~ gamma == .(s[4]))
cols <- c("#1f4e79", "#c00000", "#2e7d32", "#7b1fa2", "#ef6c00", "#455a64", "#8c564b")
y <- seq(0.0005, 0.9995, length.out = 3000)
## grid that is dense near 0 and 1 (some densities are unbounded there)
yy <- sort(unique(c(exp(seq(log(1e-12), log(0.01), length.out = 400)), seq(0.01, 0.99, length.out = 3000),
                    1 - exp(seq(log(0.01), log(1e-12), length.out = 400)))))

## ---- Figure 2: (a) pdf and (b) cdf ----
## set 5 is bimodal (first mode near y = 5e-5, shown in the inset); set 7 is increasing
pdfsets <- list(c(1.2, 1.2, 2, 1.5), c(0.5, 8, 2, 1.2), c(3, 0.2, 0.5, 2), c(0.5, 0.5, 0.5, 0.5),
                c(2.801, 0.0066, 0.0574, 1.982), c(2, 5, 1, 0.5), c(0.5, 0.5, 2, 1.6))
png("results/figure2_pdf_cdf.png", width = 3.6, height = 6.3, units = "in", res = 300)
layout(matrix(1:3, 3), heights = c(1, 1, 0.45))
par(mar = c(4, 4, 2, 1), cex = 0.7)
plot(NA, xlim = c(0, 1), ylim = c(0, 6), xlab = "y", ylab = "g(y)", main = "pdf")
for (i in seq_along(pdfsets)) { s <- pdfsets[[i]]
  lines(yy, dBOIPW(yy, s[1], s[2], s[3], s[4]), col = cols[i], lty = i, lwd = 1.3) }
mtext("(a)", side = 3, adj = -0.12, line = 0.5, font = 2)
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "y", ylab = "G(y)", main = "cdf")
for (i in seq_along(pdfsets)) { s <- pdfsets[[i]]
  lines(c(0, yy, 1), c(0, pBOIPW(yy, s[1], s[2], s[3], s[4]), 1), col = cols[i], lty = i, lwd = 1.3) }
mtext("(b)", side = 3, adj = -0.12, line = 0.5, font = 2)
par(mar = c(0, 0, 0, 0)); plot.new()
legend("center", legend = sapply(pdfsets, function(s) as.expression(lab(s))),
       col = cols, lty = 1:7, bty = "n", cex = 0.75)
## inset: the bimodal density (set 5) with y on a log scale
par(fig = c(0.36, 0.84, 0.80, 0.955), new = TRUE, mar = c(1.6, 1.8, 1.1, 0.4), cex = 0.5, mgp = c(1, 0.3, 0), tcl = -0.2)
yl <- exp(seq(log(1e-7), log(0.999), length.out = 800)); s <- pdfsets[[5]]
plot(yl, dBOIPW(yl, s[1], s[2], s[3], s[4]), type = "l", log = "x", col = cols[5], lty = 5,
     xlab = "", ylab = "", main = "set 5, y on a log scale", ylim = c(0, 2.6))
dev.off()

## ---- Figure 3: hazard rate functions ----
hazsets <- list(c(1.2, 1.2, 2, 1.5), c(3, 0.2, 0.5, 2), c(0.5, 0.5, 0.5, 0.5),
                c(1.6, 1.6, 1.4, 0.8), c(0.111, 1.291, 0.148, 2.174))
png("results/figure3_hazard.png", width = 3.6, height = 2.9, units = "in", res = 300)
par(mar = c(4, 4, 2, 1), cex = 0.7)
plot(NA, xlim = c(0, 1), ylim = c(0, 12), xlab = "y", ylab = "h(y)", main = "hazard rate function")
for (i in seq_along(hazsets)) { s <- hazsets[[i]]
  lines(y, hBOIPW(y, s[1], s[2], s[3], s[4]), col = cols[i], lty = i, lwd = 1.3) }
legend("topleft", legend = sapply(hazsets, function(s) as.expression(lab(s))),
       col = cols, lty = 1:5, bty = "n", cex = 0.7)
dev.off()

## ---- Figure 4: random numbers by the inverse transform method ----
set.seed(2024)
png("results/figure4_random.png", width = 3.6, height = 4.7, units = "in", res = 300)
par(mfrow = c(2, 1), mar = c(4, 4, 2, 1), cex = 0.7)
for (k in 1:2) {
  s <- list(c(1.2, 1.2, 2, 1.5), c(1.6, 1.6, 1.4, 0.8))[[k]]
  r <- rBOIPW(1000, s[1], s[2], s[3], s[4])
  hist(r, breaks = 30, freq = FALSE, col = "grey85", border = "grey40", xlim = c(0, 1),
       ylim = c(0, 3), main = lab(s), xlab = "y", ylab = "density")
  lines(y, dBOIPW(y, s[1], s[2], s[3], s[4]), lwd = 1.4)
  mtext(c("(a)", "(b)")[k], side = 3, adj = -0.12, line = 0.5, font = 2)
}
dev.off()
cat("Table 1 and Figures 2-4 written to results/\n")
