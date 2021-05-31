##### Create Graphs for Static Forecast Data for Manuscript

### Load in libraries and data
library(scales)
library(xtable)
setwd("ForecastApp/")
source("global.R")
load("data/forecast_fits.Rdata")
load("data/forecast_fits_isolated70.Rdata")

actual_dat <- read_csv("data/COVID-19_County_Statistics_HPSC_Ireland_(Point_Geometry)_.csv") %>% 
  filter(CountyName == "Dublin") %>%
  select(TimeStampDate, ConfirmedCovidCases)
actual_dat <- actual_dat[-(1:2), ]
actual_dat$Date <- as.Date(substr(actual_dat$TimeStampDate, 1, 10), format = "%Y/%m/%d")
actual_dat$Cases <- c(0, diff(actual_dat$ConfirmedCovidCases))
actual_dat_trim <- filter(actual_dat, (Date >= as.Date("2021-01-26")) & (Date <= as.Date("2021-03-28")))

### Sum across compartments
comp_sel <- function(x, y) {
  grepl(x, y)
}
comps <- c("S_", "Ev_", "Ip_", "IA_", "Ii_", "It_", "Iti_", "Iq_", "R_")
comp_inds <- lapply(comps, function(x, y) grepl(x, y), y = colnames(forecast_fits[[1]][[1]]))

sum_comp_inds <- function(x, y) {
  rowSums(y[x])
}

### Define Scenarios
scenarios <- c("No Intervention -> No Intervention",
               "Level 3 -> Level 3",
               "Level 5 -> Level 5",
               "Level 3 -> Level 5",
               "Level 2 -> Level 5",
               "Level 1 -> Level 5")

### Calculate deaths
age_deaths <- function(x, forecaster = FALSE, just_cases = FALSE, back = 63) {
  denom_1 <- def_pars["Dv"] - def_pars["Cv"] + def_pars["L"]
  denom_2 <- denom_1 - def_pars["TT"]
  N <- nrow(x)
  get_rows <- (N-back):N
  moved_removed <- (x[get_rows, 66:81] + x[get_rows, 114:129])/denom_1 + x[get_rows, 98:113]/denom_2
  if (!forecaster) {
    ags <- matrix(0, nrow = nrow(moved_removed), ncol = 8)
    ags[, 1] <- rowSums(moved_removed[, 1:3])
    ags[, 2] <- rowSums(moved_removed[, 4:5])
    ags[, 3] <- rowSums(moved_removed[, 6:7])
    ags[, 4] <- rowSums(moved_removed[, 8:9])
    ags[, 5] <- rowSums(moved_removed[, 10:11])
    ags[, 6] <- rowSums(moved_removed[, 12:13])
    ags[, 7] <- rowSums(moved_removed[, 14:15])
    ags[, 8] <- moved_removed[, 16]
    if ( just_cases ) { return( ags ) }
    return(ags %*% est_deaths$Estimated_Deaths)
  }
  ags <- matrix(0, nrow = nrow(moved_removed), ncol = 9)
  ags[, 1] <- rowSums(moved_removed[, 1:3])
  ags[, 2] <- rowSums(moved_removed[, 4:5])
  ags[, 3] <- rowSums(moved_removed[, 6:7])
  ags[, 4] <- rowSums(moved_removed[, 8:9])
  ags[, 5] <- rowSums(moved_removed[, 10:11])
  ags[, 6] <- rowSums(moved_removed[, 12:13])
  ags[, 7] <- moved_removed[, 14]
  ags[, 8] <- moved_removed[, 15]
  ags[, 9] <- moved_removed[, 16]
  est_alt <- est_deaths$Estimated_Deaths[c(1:8, 8)]
  prop_70 <- as.numeric( dub_population[15, 2] / sum(dub_population[14:15, 2]) )
  est_alt[7] <- est_alt[7] * (1 - prop_70)
  est_alt[8] <- est_alt[8] * prop_70
  if ( just_cases ) { return( ags ) }
  ags %*% est_alt
}

make_graphs <- function(dat, filename, forecaster = FALSE) {
  
  len <- length(dat[[1]])
  UL <- MID <- LL <- vector("list", len)
  for (i in 1:len) {
    UL[[i]] <- as.data.frame( sapply(comp_inds, sum_comp_inds, y = dat$UL[[i]]) )
    MID[[i]] <- as.data.frame( sapply(comp_inds, sum_comp_inds, y = dat$MID[[i]]) )
    LL[[i]] <- as.data.frame( sapply(comp_inds, sum_comp_inds, y = dat$LL[[i]]) )
    colnames(UL[[i]]) <- colnames(MID[[i]]) <- colnames(LL[[i]]) <- comps
  }
  
  ### Add in 'All', 'Symptomatic' compartments
  add_comps <- function(x) {
    x$Sy <- x$Ii_ + x$It_ + x$Iti_ + x$Iq_
    x$Al <- x$Ip_ + x$IA + x$Sy
    x
  }
  
  UL <- lapply(UL, add_comps)
  MID <- lapply(MID, add_comps)
  LL <- lapply(LL, add_comps)
    
  ### Draw Plots
  start_date <- as.Date("2021-02-01")
  xax <- seq.Date(start_date - 7, start_date + 7*8, 7)
  xval <- seq.Date(min(xax), max(xax), 1)
  xval_disp <- as.Date(c("2021-01-28", "2021-03-26"))
  n_date <- as.numeric( difftime(max(xax), min(xax), units = "days") )
  
  # Draw empty plot that we will add to
  truncate_dat <- function(x) {
    x[["It_"]][(nrow(x) - n_date):nrow(x)] / def_pars["TT"]
  }
  
  UL_all <- lapply(UL, truncate_dat)
  MID_all <- lapply(MID, truncate_dat)
  LL_all <- lapply(LL, truncate_dat)
  
  UL_y <- UL_all[-c(1, 7)]
  MID_y <- MID_all[-c(1, 7)]
  LL_y <- LL_all[-c(1, 7)]
  max_ind <- which.max( sapply(UL_y, max) )
  
  pdf(filename)
  plot(x = xval, y = UL_y[[max_ind]], type = "n", ylab = "Expected Daily Cases", xaxt = "n",
       xlab = "", ylim = c(0, 6500), 
       xlim = xval_disp,
       main = "")
  axis.Date(1, at = xax, las = 2)
  grid()
  
  # Add forecasts
  use_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  N <- length(MID_y)
  for (i in 1:N) {
    polygon(c(xval,rev(xval)),c(LL_y[[i]],rev(UL_y[[i]])),col=alpha(use_cols[i], 0.1),
          border = alpha(use_cols[i], 0.1))
  }
  
  for (i in 1:N) {
    lines(x = xval, y = MID_y[[i]], col = use_cols[i])
  }
  
  lines(Cases ~ Date, data = actual_dat_trim, col = "#999999")
  
  legend("topleft", col = c(use_cols[1:N], "#999999"), lty = 1, 
         legend = c(scenarios[2:6], "Actual Cases"),
         bty = "n")
  dev.off()
  
  UL_deaths <- Map(age_deaths, forecast_fits$UL, forecaster = forecaster)
  MID_deaths <- Map(age_deaths, forecast_fits$MID, forecaster = forecaster)
  LL_deaths <- Map(age_deaths, forecast_fits$LL, forecaster = forecaster)
  
  UL_y <- UL_deaths[-c(1, 7)]
  MID_y <- MID_deaths[-c(1, 7)]
  LL_y <- LL_deaths[-c(1, 7)]
  max_ind <- which.max( sapply(UL_y, max) )
  
  pdf(paste0("deaths_", filename))
  plot(x = xval, y = UL_y[[max_ind]], type = "n", ylab = "Expected Daily Deaths", xaxt = "n",
       xlab = "", ylim = c(0, 50), 
       xlim = xval_disp,
       main = "")
  axis.Date(1, at = xax, las = 2)
  grid()
  
  # Add forecasts
  use_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7", "#999999")
  
  N <- length(MID_y)
  for (i in 1:N) {
    polygon(c(xval,rev(xval)),c(LL_y[[i]],rev(UL_y[[i]])),col=alpha(use_cols[i], 0.1),
          border = alpha(use_cols[i], 0.1))
  }
  
  for (i in 1:N) {
    lines(x = xval, y = MID_y[[i]], col = use_cols[i])
  }
  
  legend("topleft", col = use_cols[1:N], lty = 1, 
         legend = scenarios[2:6],
         bty = "n")
  dev.off()
  
}

make_graphs(forecast_fits, "forecast_fits.pdf")
make_graphs(forecast_fits_isolated70, "forecast_fits_isolated70.pdf", forecaster = TRUE)
  
### Projection Plots
policy <- c('No Intervention', 'Level 3', 'Level 5', 'Level 3 > Level 5',
            'Level 2 > Level 5', 'Level 1 > Level 5')
UL_deaths <- round( sapply( Map(age_deaths, forecast_fits$UL, back = 56), sum ), 0 )[-7]
MID_deaths <- round( sapply( Map(age_deaths, forecast_fits$MID, back = 56), sum ), 0 )[-7]
LL_deaths <- round( sapply( Map(age_deaths, forecast_fits$LL, back = 56), sum ), 0 )[-7]

UL_cases <- round( sapply( Map(age_deaths, forecast_fits$UL, just_cases = TRUE, back = 56), sum ), 0 )[-7]
MID_cases <- round( sapply( Map(age_deaths, forecast_fits$MID, just_cases = TRUE, back = 56), sum ), 0 )[-7]
LL_cases <- round( sapply( Map(age_deaths, forecast_fits$LL, just_cases = TRUE, back = 56), sum ), 0 )[-7]

projection_df <- data.frame(Policy = policy, 
                            Cases = paste0(MID_cases, ' (', LL_cases, ', ', UL_cases, ')'),
                            Deaths = paste0(MID_deaths, ' (', LL_deaths, ', ', UL_deaths, ')'))

print(xtable(projection_df), include.rownames = FALSE)

### Projection Plots (with isolation)
UL_deaths2 <- round( sapply( Map(age_deaths, forecast_fits_isolated70$UL, forecaster = TRUE, back = 56), sum ), 0 )[-7]
MID_deaths2 <- round( sapply( Map(age_deaths, forecast_fits_isolated70$MID, forecaster = TRUE, back = 56), sum ), 0 )[-7]
LL_deaths2 <- round( sapply( Map(age_deaths, forecast_fits_isolated70$LL, forecaster = TRUE, back = 56), sum ), 0 )[-7]

UL_cases2 <- round( sapply( Map(age_deaths, forecast_fits_isolated70$UL, forecaster = TRUE, just_cases = TRUE, back = 56), sum ), 0 )[-7]
MID_cases2 <- round( sapply( Map(age_deaths, forecast_fits_isolated70$MID, forecaster = TRUE, just_cases = TRUE, back = 56), sum ), 0 )[-7]
LL_cases2 <- round( sapply( Map(age_deaths, forecast_fits_isolated70$LL, forecaster = TRUE, just_cases = TRUE, back = 56), sum ), 0 )[-7]

projection_df2 <- data.frame(Policy = policy, 
                             Cases = paste0(MID_cases2, ' (', LL_cases2, ', ', UL_cases2, ')'),
                             Deaths = paste0(MID_deaths2, ' (', LL_deaths2, ', ', UL_deaths2, ')'))

print(xtable(projection_df2), include.rownames = FALSE)

### Plots
all_cases <- sum( actual_dat_trim$Cases[7:nrow(actual_dat_trim)] )
pdf('case_projections.pdf')
par(mar = c(8, 4, 4, 2) + 0.1)
plot(MID_cases[-1], ylim = c(0, 200000), col = 'white',
     ylab = 'Expected Cases', xlab = '', xaxt = 'n', yaxt = 'n')
axis(1, at = 1:5, labels = policy[-1], las = 2)
axis(2, at = seq(0, 200000, 50000), labels = seq(0, 20, 5))
mtext(expression(x10^4), at = c(1, 20000))
grid()
arrows(x0 = 1:5, y0 = LL_cases[-1], y1 = UL_cases[-1], angle = 90,
       code = 3, length = 0.1, col = 'dodgerblue')
points(MID_cases[-1], pch = '-', cex = 2, col = 'dodgerblue')
abline(h = all_cases, lty = 2, col = 'orange')
dev.off()

pdf('case_projections_isolated.pdf')
par(mar = c(8, 4, 4, 2) + 0.1)
plot(MID_cases2[-1], ylim = c(0, 200000), col = 'white',
     ylab = 'Expected Cases', xlab = '', xaxt = 'n', yaxt = 'n')
axis(1, at = 1:5, labels = policy[-1], las = 2)
axis(2, at = seq(0, 200000, 50000), labels = seq(0, 20, 5))
mtext(expression(x10^4), at = c(1, 20000))
grid()
arrows(x0 = 1:5, y0 = LL_cases2[-1], y1 = UL_cases2[-1], angle = 90,
       code = 3, length = 0.1, col = 'dodgerblue')
points(MID_cases2[-1], pch = '-', cex = 2, col = 'dodgerblue')
abline(h = all_cases, lty = 2, col = 'orange')
dev.off()

### Actual Case Plot
pdf("Case_Plot.pdf")
plot(Cases ~ Date, data = actual_dat_trim, type = "l", ylim = c(0, max(Cases)),
     col = 'dodgerblue')
grid()
dev.off()
