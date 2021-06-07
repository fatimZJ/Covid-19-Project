##### Create Graphs for Static Forecast Data for Manuscript

### Load in libraries and data
library(scales)
library(xtable)
setwd("ForecastApp/")
source("global.R")
load("data/forecast_fits.Rdata")
load("data/forecast_fits_isolated70.Rdata")

start_date <- as.Date("2021-02-01")
start_date2 <- start_date + (4*7) + 1
date_start_seq <- c(start_date, start_date2)
date_start_end <- c(start_date2 - 1, start_date2 + (4*7))

date_end <- date_start_end[2]

actual_dat <- read_csv("data/COVID-19_County_Statistics_HPSC_Ireland_(Point_Geometry)_.csv") %>% 
  filter(CountyName == "Dublin") %>%
  select(TimeStampDate, ConfirmedCovidCases)
actual_dat <- actual_dat[-(1:2), ]
actual_dat$Date <- as.Date(substr(actual_dat$TimeStampDate, 1, 10), format = "%Y/%m/%d")
actual_dat$Cases <- c(1, diff(actual_dat$ConfirmedCovidCases))
actual_dat_trim <- filter(actual_dat, (Date >= as.Date("2021-01-26")) & (Date <= date_end))

### Sum across compartments
comp_sel <- function(x, y) {
  grepl(x, y)
}
comps <- c("S_", "Ev_", "Ip_", "IA_", "Ii_", "It_", "Iti_", "Iq_", "R_", "Cc_")
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
  
  N <- nrow(x)
  get_rows <- (N-back):N
  if (just_cases) {
    return ( rowSums(x[get_rows, grepl("It_", names(x))] / def_pars["TT"]) )
  }
  
  denom_1 <- def_pars["Dv"] - def_pars["Cv"] + def_pars["L"]
  denom_2 <- denom_1 - def_pars["TT"]
  
  Removed <- (x[get_rows, grepl("Ii_", names(x))] + x[get_rows, grepl("Iq_", names(x))])/denom_1 + 
    x[get_rows, grepl("Iti_", names(x))]/denom_2
  
  if (!forecaster) {
    ags <- matrix(0, nrow = nrow(Removed), ncol = 8)
    ags[, 1] <- rowSums(Removed[, 1:3])
    ags[, 2] <- rowSums(Removed[, 4:5])
    ags[, 3] <- rowSums(Removed[, 6:7])
    ags[, 4] <- rowSums(Removed[, 8:9])
    ags[, 5] <- rowSums(Removed[, 10:11])
    ags[, 6] <- rowSums(Removed[, 12:13])
    ags[, 7] <- rowSums(Removed[, 14:15])
    ags[, 8] <- Removed[, 16]
    return(ags %*% est_deaths$Estimated_Deaths)
  }
  
  ags <- matrix(0, nrow = nrow(Removed), ncol = 9)
  ags[, 1] <- rowSums(Removed[, 1:3])
  ags[, 2] <- rowSums(Removed[, 4:5])
  ags[, 3] <- rowSums(Removed[, 6:7])
  ags[, 4] <- rowSums(Removed[, 8:9])
  ags[, 5] <- rowSums(Removed[, 10:11])
  ags[, 6] <- rowSums(Removed[, 12:13])
  ags[, 7] <- Removed[, 14]
  ags[, 8] <- Removed[, 15]
  ags[, 9] <- Removed[, 16]
  
  est_alt <- est_deaths$Estimated_Deaths[c(1:8, 8)]
  prop_70 <- as.numeric( dub_population[15, 2] / sum(dub_population[14:15, 2]) )
  est_alt[7] <- est_alt[7] * (1 - prop_70)
  est_alt[8] <- est_alt[8] * prop_70
  
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
  xval_disp <- as.Date(c("2021-01-28", "2021-03-27"))
  n_date <- as.numeric( difftime(max(xax), min(xax), units = "days") )
  
  # Draw empty plot that we will add to
  truncate_dat <- function(x) {
    x[["It_"]][(nrow(x) - n_date):nrow(x)] / def_pars["TT"]
  }
  
  UL_all <- lapply(UL, truncate_dat)
  MID_all <- lapply(MID, truncate_dat)
  LL_all <- lapply(LL, truncate_dat)
  
  UL_y <- UL_all[2:5]
  MID_y <- MID_all[2:5]
  LL_y <- LL_all[2:5]
  max_ind <- which.max( sapply(UL_y, max) )
  
  pdf(filename)
  plot(x = xval, y = UL_y[[max_ind]], type = "n", ylab = "Expected Daily Cases", xaxt = "n",
       xlab = "", ylim = c(0, 900), 
       xlim = xval_disp,
       main = "")
  axis.Date(1, at = xax, las = 2)
  grid()
  
  # Add forecasts
  use_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#999999")
  
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
         legend = c(scenarios[2:5], "Actual Cases"),
         bty = "n")
  dev.off()
  
  UL_deaths <- Map(age_deaths, forecast_fits$UL, forecaster = forecaster)
  MID_deaths <- Map(age_deaths, forecast_fits$MID, forecaster = forecaster)
  LL_deaths <- Map(age_deaths, forecast_fits$LL, forecaster = forecaster)
  
  UL_y <- UL_deaths[2:5]
  MID_y <- MID_deaths[2:5]
  LL_y <- LL_deaths[2:5]
  max_ind <- which.max( sapply(UL_y, max) )
  
  pdf(paste0("deaths_", filename))
  plot(x = xval, y = UL_y[[max_ind]], type = "n", ylab = "Expected Daily Deaths", xaxt = "n",
       xlab = "", ylim = c(0, 3.5), 
       xlim = xval_disp,
       main = "")
  axis.Date(1, at = xax, las = 2)
  grid()
  
  # Add forecasts
  use_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#999999")
  
  N <- length(MID_y)
  for (i in 1:N) {
    polygon(c(xval,rev(xval)),c(LL_y[[i]],rev(UL_y[[i]])),col=alpha(use_cols[i], 0.1),
          border = alpha(use_cols[i], 0.1))
  }
  
  for (i in 1:N) {
    lines(x = xval, y = MID_y[[i]], col = use_cols[i])
  }
  
  legend("topleft", col = use_cols[1:N], lty = 1, 
         legend = scenarios[2:5],
         bty = "n")
  dev.off()
  
}

make_graphs(forecast_fits, "forecast_fits.pdf")
make_graphs(forecast_fits_isolated70, "forecast_fits_isolated70.pdf", forecaster = TRUE)
  
### Projection Plots
policy <- c('No Intervention', 'Level 3', 'Level 5', 'Level 3 > Level 5',
            'Level 2 > Level 5')
go_back <- (difftime(date_end, start_date)-1)[[1]]
UL_deaths <- sapply( Map(age_deaths, forecast_fits$UL, back = go_back), sum )[1:5]
MID_deaths <- sapply( Map(age_deaths, forecast_fits$MID, back = go_back), sum )[1:5]
LL_deaths <- sapply( Map(age_deaths, forecast_fits$LL, back = go_back), sum )[1:5]

UL_cases <- sapply( Map(age_deaths, forecast_fits$UL, just_cases = TRUE, back = go_back), sum )[1:5]
MID_cases <- sapply( Map(age_deaths, forecast_fits$MID, just_cases = TRUE, back = go_back), sum )[1:5]
LL_cases <- sapply( Map(age_deaths, forecast_fits$LL, just_cases = TRUE, back = go_back), sum )[1:5]

projection_df <- data.frame(Policy = policy, 
                            Cases = paste0(round(MID_cases/100, 0), ' (', 
                                           round(LL_cases/100, 0), ', ', 
                                           round(UL_cases/100, 0), ')'),
                            Deaths = paste0(round(MID_deaths, 0), ' (', 
                                            round(LL_deaths, 0), ', ', 
                                            round(UL_deaths, 0), ')'))

print(xtable(projection_df), include.rownames = FALSE)

### Projection Plots (with isolation)
UL_deaths2 <- sapply( Map(age_deaths, forecast_fits_isolated70$UL, forecaster = TRUE, back = go_back), sum )[1:5]
MID_deaths2 <- sapply( Map(age_deaths, forecast_fits_isolated70$MID, forecaster = TRUE, back = go_back), sum )[1:5]
LL_deaths2 <- sapply( Map(age_deaths, forecast_fits_isolated70$LL, forecaster = TRUE, back = go_back), sum )[1:5]

UL_cases2 <- sapply( Map(age_deaths, forecast_fits_isolated70$UL, forecaster = TRUE, just_cases = TRUE, back = go_back), sum )[1:5]
MID_cases2 <- sapply( Map(age_deaths, forecast_fits_isolated70$MID, forecaster = TRUE, just_cases = TRUE, back = go_back), sum )[1:5]
LL_cases2 <- sapply( Map(age_deaths, forecast_fits_isolated70$LL, forecaster = TRUE, just_cases = TRUE, back = go_back), sum )[1:5]

projection_df2 <- data.frame(Policy = policy, 
                             Cases = paste0(round(MID_cases2/100, 0), ' (', 
                                            round(LL_cases2/100, 0), ', ', 
                                            round(UL_cases2/100, ), ')'),
                             Deaths = paste0(round(MID_deaths2, 0), ' (', 
                                             round(LL_deaths2, 0), ', ', 
                                             round(UL_deaths2, 0), ')'))

print(xtable(projection_df2), include.rownames = FALSE)

### Plots
ind_start <- which(actual_dat_trim$Date == start_date)
ind_end <- which(actual_dat_trim$Date == date_end)
all_cases <- sum( actual_dat_trim$Cases[ind_start:ind_end] )
actual_dat_trim$Date[ind_start:ind_end]

pdf('case_projections.pdf')
par(mar = c(8, 4, 4, 2) + 0.1)
plot(MID_cases[-1], ylim = c(0, 35000), col = 'white',
     ylab = 'Expected Cases', xlab = '', xaxt = 'n')#, yaxt = 'n')
axis(1, at = 1:4, labels = policy[-1], las = 2)
#axis(2, at = seq(0, 300000, 50000), labels = seq(0, 30, 5))
#mtext(expression(x10^4), at = c(1, 20000))
grid()
arrows(x0 = 1:4, y0 = LL_cases[-1], y1 = UL_cases[-1], angle = 90,
       code = 3, length = 0.1, col = 'dodgerblue')
points(MID_cases[-1], pch = 20, cex = 1, col = 'dodgerblue')
abline(h = all_cases, lty = 2, col = 'orange')
dev.off()

pdf('case_projections_isolated.pdf')
par(mar = c(8, 4, 4, 2) + 0.1)
plot(MID_cases2[-1], ylim = c(0, 35000), col = 'white',
     ylab = 'Expected Cases', xlab = '', xaxt = 'n')#, yaxt = 'n')
axis(1, at = 1:4, labels = policy[-1], las = 2)
#axis(2, at = seq(0, 300000, 50000), labels = seq(0, 30, 5))
#mtext(expression(x10^4), at = c(1, 20000))
grid()
arrows(x0 = 1:4, y0 = LL_cases2[-1], y1 = UL_cases2[-1], 
       angle = 90, code = 3, length = 0.1, col = 'dodgerblue')
points(MID_cases2[-1], pch = 20, cex = 1, col = 'dodgerblue')
abline(h = all_cases, lty = 2, col = 'orange')
dev.off()

### Actual Case Plot
pdf("Case_Plot.pdf")
plot(Cases ~ Date, data = actual_dat_trim, type = "l", ylim = c(0, max(Cases)),
     col = 'dodgerblue')
grid()
dev.off()
