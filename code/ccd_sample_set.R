##### CCD design specs (using package by Jason)
library(ccdpts)

# Build intervals matrix using 2 standard deviations
parnames <- c("tau_C", "tau_P", "tau_D_SC",
              "tau_D_C", "R0", "alpha", "fa")
means <- c(1.63, 0.69, 1.63, 2.57, 1.3, 0.55, -1.6)
sds <- c(0.5, 0.35, 0.53, 0.09, 0.15, 0.15, 0.1)

inter <- matrix(c(means - 2*sds, means, means + 2*sds), ncol = 3,
                dimnames = list(parnames, c("Lower", "Mean", "Upper")))

# CCD design
ccd_set <- ccd.design.grid(intervals = inter, range = 2)
set_mat <- data.frame(t(ccd_set$points))

fun_list <- list(exp, exp, exp, exp, exp, identity, exp)
trans_set_mat <- mapply(function(x, y) {x(y)}, x = fun_list, y = set_mat)
colnames(trans_set_mat) <- parnames

ccd_set$weights

