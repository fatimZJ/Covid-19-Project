##### CCD design specs (using package by Jason)
library(ccdpts)

# Build intervals matrix using 2 standard deviations
parnames <- c("tau_C", "tau_P", "tau_D_SC",
              "tau_D_C", "R0", "alpha", "k",
              "fa", "ft", "fk", "T")
means <- c(1.63, 0.69, 1.63, 2.57, 1.3, 0.55, 0.05, -1.6, 0.75, 0.05, 2)
sds <- c(0.5, 0.35, 0.53, 0.09, 0.15, 0.15, 0.025, 0.1, 0.075, 0.025, 0.5)
inter <- matrix(c(means - 2*sds, means, means + 2*sds), ncol = 3,
                dimnames = list(parnames, c("Lower", "Mean", "Upper")))

# CCD design
ccd_set <- ccd.design.grid(intervals = inter , range = 2)
set_mat <- ccd_set$points
ccd_set$weights

