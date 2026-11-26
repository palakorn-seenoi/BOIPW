# BOIPW: the bounded odd inverse Pareto Weibull distribution in R

R code accompanying the paper

> Seenoi P, Saengthong P. Bounded odd inverse Pareto Weibull distribution: properties, estimation, and applications. *Journal of Applied Research on Science and Technology (JARST)*. Submitted, 2026.

The code uses **base R only** (tested with R 4.3.3); no extra packages are needed.

Repository: <https://github.com/palakorn-seenoi/BOIPW>

```sh
git clone https://github.com/palakorn-seenoi/BOIPW.git
cd BOIPW
```

Or download the ZIP file from the green **Code** button on GitHub and unzip it.

## The distribution

For $0<y<1$ and $\alpha,\beta,\theta,\gamma>0$, let $K(y)=\exp\{-\theta(-\ln y)^{\gamma}\}$ and $\bar K(y)=1-K(y)$. The BOIPW cdf is

$$G(y)=1-\left[\frac{\bar K(y)}{\beta+(1-\beta)\bar K(y)}\right]^{\alpha}.$$

Special cases: BOIPE ($\gamma=1$), BMOEW ($\alpha=1$), UEW, unit exponentiated Weibull ($\beta=1$), BMOEE ($\alpha=\gamma=1$), Kumaraswamy ($\beta=\gamma=1$) and unit-Weibull ($\alpha=\beta=1$).

## Repository structure

```
R/boipw.R                  d/p/q/r/h functions, moments, sub-models
R/fit.R                    maximum likelihood, K-S / A-D / Cramer-von Mises tests, LR tests
data/ncyoungmale.csv       Data I  (n = 90), see data/README.md
data/catholic.csv          Data II (n = 46), see data/README.md
scripts/00_verify_formulas.R  numerical checks of all analytical results
scripts/01_properties.R    Table 1, Figures 2-4
scripts/02_simulation.R    Table 2 and Tables S1-S2 of the Supporting Information (Monte Carlo study)
scripts/03_applications.R  Tables 3-10, Figures 5-6
run_all.R                  runs all four scripts
results/                   output produced by the scripts
Supporting_Information.pdf Supporting Information of the paper (Tables S1-S2, additional simulation results)
Supporting_Information.docx the same file in Word format
```

## Quick start

```r
source("R/boipw.R"); source("R/fit.R")

dBOIPW(0.5, alpha = 1.2, beta = 1.2, theta = 2, gamma = 1.5)   # pdf
pBOIPW(0.5, 1.2, 1.2, 2, 1.5)                                  # cdf
qBOIPW(c(0.25, 0.5, 0.75), 1.2, 1.2, 2, 1.5)                   # quantiles
hBOIPW(0.5, 1.2, 1.2, 2, 1.5)                                  # hazard rate
y <- rBOIPW(1000, 1.2, 1.2, 2, 1.5)                            # random numbers
summaryBOIPW(1.2, 1.2, 2, 1.5)                                 # mean, median, variance, skewness, kurtosis

# fit BOIPW and a sub-model to data, goodness of fit and LR test
y <- read.csv("data/ncyoungmale.csv")$proportion
f4 <- fit_model("BOIPW", y)
f3 <- fit_model("BOIPE", y)
f4$par; f4$se; f4$AIC; f4$BIC
gof_stats(f4, y)         # K-S, A-D and Cramer-von Mises statistics with p-values
lr_test(f4, f3)          # likelihood ratio test of gamma = 1
```

Available models in `fit_model()`: the BOIPW distribution and its sub-models `"BOIPW"`, `"BOIPE"`, `"BMOEW"`, `"UEW"` (unit exponentiated Weibull, beta = 1), `"BMOEE"`, `"Kumaraswamy"` and `"UW"` (unit-Weibull), and the non-nested competitors `"EKw"` (exponentiated Kumaraswamy), `"UBXII"` (unit Burr XII) and `"Beta"`. Each log-parameter is restricted to [-12, 12]; a fit on the edge of this region is flagged with `at_bound = TRUE` and its standard errors are `NA`.

## Reproducing the paper

From the repository root:

```sh
Rscript scripts/00_verify_formulas.R # about 2 seconds
Rscript scripts/01_properties.R      # about 5 seconds
Rscript scripts/03_applications.R    # about 10 seconds
Rscript scripts/02_simulation.R 1000 4   # 1000 replications on 4 cores, about 10-15 minutes
```

On Windows use one core (`Rscript scripts/02_simulation.R 1000 1`), which takes about four times as long. In the simulation, each sample is fitted from 13 starting points (the true values, two perturbations of them, (1, 1, 1, 1) and nine random points), because the BOIPW likelihood is often bimodal; `results/tableS1_S2_simulation.csv` (Tables S1 and S2 of the Supporting Information) reports the percentage of boundary estimates, the coverage of the 95% Wald intervals, how often the three starts based on the true values (the true vector and two perturbations of it) would have missed the global maximum, how often the ten starts that do not use the true values reach the same maximum, the coverage of Wald intervals on the log scale, and the mean Kolmogorov distance between the fitted and the true cdfs. The simulation restricts each log-parameter to [-10, 10] and the applications to [-12, 12]; interior estimates are the same under both. Figure 1 is a diagram and is not produced by the code.

### Numerical verification of the formulas

`scripts/00_verify_formulas.R` compares every closed-form result of the paper with an independent numerical computation and prints the largest error of each check (also saved to `results/verification.csv`):

* pdf (5) versus the numerical derivative of the cdf (4), hazard rate (7), and total probability 1 (Proposition 1, Remark 1);
* the boundary behaviour of the pdf and hazard rate (Proposition 2);
* the quantile function (8), through G(Q(p)) = p (Proposition 3);
* the three moment expressions (9) (Proposition 4);
* the geometric and negative binomial expansions (10)-(11) and the moment series (12)-(13) (Proposition 5);
* the stress-strength reliability (16), by numerical integration and by simulation (Proposition 6);
* the score functions (20)-(23) and the second derivatives (B1)-(B8), against finite differences on both data sets.

The script stops with an error if any check exceeds its tolerance.

### Agreement with the paper

* **Table 1 and Tables 3-10** (real data): the R results agree with the paper to the reported number of decimals, including the log-likelihoods, AIC, BIC, the K-S, A-D and Cramer-von Mises statistics and the LR statistics. The standard errors are computed from the analytic second derivatives of the log-likelihood (Appendix B of the paper; functions `boipw_score()`, `boipw_hessian()` and `loglik_hessian()` in `R/fit.R`), so they do not depend on numerical differentiation.
* **Table 2 and Tables S1-S2** (simulation): the results agree up to Monte Carlo error but not to the last digit, because the random numbers differ from those of the original run. The mean bias and RMSE for small $n$ are dominated by a few extreme estimates (see the paper), so they vary more between runs than the median bias.

### Notes on the tests

* K-S: exact two-sided p-value (`ks.test(..., exact = TRUE)`).
* A-D: p-value from Marsaglia and Marsaglia (2004).
* Cramer-von Mises: p-value from the finite-sample approximation of Csörgő and Faraway (1996).

## License

MIT, see [LICENSE](LICENSE).
