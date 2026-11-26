## Reproduce every table and figure of the paper (except the diagram in Figure 1).
## Run from the repository root:  Rscript run_all.R
## The simulation (1000 replications) takes a while on one core; on Linux/macOS
## you can run scripts/02_simulation.R separately with more cores (see README.md).
source("scripts/00_verify_formulas.R")
source("scripts/01_properties.R")
source("scripts/03_applications.R")
commandArgs <- function(...) c("1000", "1")
source("scripts/02_simulation.R")
