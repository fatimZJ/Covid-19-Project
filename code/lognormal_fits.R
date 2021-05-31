xx <- seq(0.01, 20, 0.01)
yy <- dlnorm(xx, 2, 1)
plot(xx, yy)

### Pre-symptmatic period
fit.lognormal <- function(pars, quantiles = c(2.5, 97.5), pos = c(1, 4)) {
  
  yy <- qlnorm(quantiles/100, pars[1], pars[2])
  sum( (yy - pos)^2 )
  
}

optim(par = c(0, 1), fit.lognormal)

### Asymptomatic Infectious Period
xx <- seq(0.01, 0.99, 0.01)
yy <- qgamma(xx, 3, 0.5)
plot(xx, yy)

gam_quants <- 1:99
gam_pos <- qgamma(gam_quants/100, 3, 0.5)

optim(par = c(0, 1), fit.lognormal, quantiles = gam_quants, pos = gam_pos)

yy <- qlnorm(xx, 1.66, 0.53)
points(xx, yy, col = "red")

### Symptomatic infectious period
C_quants <- c(2.5, 97.5)
C_pos <- c(10.9, 15.8)

C_opt <- optim(par = c(0, 1), fit.lognormal, quantiles = C_quants, pos = C_pos)

yy <- qlnorm(xx, C_opt$par[1], C_opt$par[2])

plot(xx, yy)

### Basic R
R_quants <- c(2.5, 97.5)
R_pos <- c(2.7664, 4.9079)

R_opt <- optim(par = c(0, 1), fit.lognormal, quantiles = R_quants, pos = R_pos)

yy <- qlnorm(xx, R_opt$par[1], R_opt$par[2])

plot(xx, yy)

mean(yy)

### Asymptomatic Proportion
fa_quants <- c(2.5, 97.5)
fa_pos <- c(0.17, 0.25)

fa_opt <- optim(par = c(0, 1), fit.lognormal, quantiles = fa_quants, pos = fa_pos)

yy <- qlnorm(xx, fa_opt$par[1], fa_opt$par[2])

plot(xx, yy)

mean(yy)

