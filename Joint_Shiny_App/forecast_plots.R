##### Create Graphs for Static Forecast Data for Manuscript

### Load in libraries and data
library(scales)
setwd("Joint_Shiny_App/")
source("global.R")
load("data/forecast_fits.Rdata")
load("data/forecast_fits_isolated70.Rdata")

actual_dat <- read_csv("data/Daily_County_Cases.csv") %>% 
  filter(CountyName == "Dublin") %>%
  select(TimeStamp, ConfirmedCovidCases)
actual_dat <- actual_dat[-1, ]
actual_dat$Date <- as.Date(substr(actual_dat$TimeStamp, 1, 10), format = "%Y/%m/%d")
actual_dat$Cases <- c(0, diff(actual_dat$ConfirmedCovidCases))
actual_dat_trim <- filter(actual_dat, (Date >= as.Date("2020-11-25")) & (Date <= as.Date("2021-01-26")))

### Sum across compartments
comp_sel <- function(x, y) {
  grepl(x, y)
}
comps <- c("S_", "Ev_", "Ip_", "IA_", "Ii_", "It_", "Iti_", "Iq_", "R_")
comp_inds <- lapply(comps, function(x, y) grepl(x, y), y = colnames(forecast_fits[[1]][[1]]))

sum_comp_inds <- function(x, y) {
  rowSums(y[x])
}

make_graphs <- function(dat, filename) {
  UL <- MID <- LL <- vector("list", 5)
  for (i in 1:5) {
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
  # Prepare x axis
  start_date <- as.Date("2020-12-01")
  xax <- seq.Date(start_date - 7, start_date + 7*8, 7)
  xval <- seq.Date(min(xax), max(xax), 1)
  xval_disp <- as.Date(c("2020-11-27", "2021-01-23"))
  n_date <- as.numeric( difftime(max(xax), min(xax), units = "days") )
  
  # Draw empty plot that we will add to
  truncate_dat <- function(x) {
    x[["It_"]][(nrow(x) - n_date):nrow(x)] / 7
  }
  
  UL_all <- lapply(UL, truncate_dat)
  MID_all <- lapply(MID, truncate_dat)
  LL_all <- lapply(LL, truncate_dat)
  
  UL_y <- UL_all[-1]
  MID_y <- MID_all[-1]
  LL_y <- LL_all[-1]
  max_ind <- which.max( sapply(UL_y, max) )
  
  pdf(filename)
  plot(x = xval, y = UL_y[[max_ind]], type = "n", ylab = "", xaxt = "n",
       xlab = "Date", ylim = c(0, max(c(UL_y[[max_ind]], actual_dat_trim$Cases))), xlim = xval_disp,
       main = "Case Forecast")
  axis.Date(1, at = xax, las = 2)
  abline(v = xax, lty = 2, col = "lightgray")
  abline(h = seq(0, 5000, 1000), lty = 2, col = "lightgray")
  
  # Add forecasts
  use_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")#, "#0072B2")
  
  N <- length(MID_y)
  for (i in 1:N) {
    polygon(c(xval,rev(xval)),c(LL_y[[i]],rev(UL_y[[i]])),col=alpha(use_cols[i], 0.1),
          border = alpha(use_cols[i], 0.1))
  }
  
  for (i in 1:N) {
    lines(x = xval, y = MID_y[[i]], col = use_cols[i])
  }
  
  lines(Cases ~ Date, data = actual_dat_trim)
  
  legend("topleft", col = c(use_cols, "#000000"), lty = 1, 
         legend = c(paste("Scenario", 2:5), "Actual Cases"),
         bty = "n")
  dev.off()
  
  # Plot no Intervention
  pdf(paste0("No_Intervention_", filename))
  plot(x = xval, y = UL_all[[1]], ylab = "", xaxt = "n",
       xlab = "Date", xlim = xval_disp, type = "n",
       ylim = c(0, 40000),
       main = "Case Forecast", col = use_cols[1])
  axis.Date(1, at = xax, las = 2)
  abline(v = xax, lty = 2, col = "lightgray")
  abline(h = seq(0, 40000, 10000), lty = 2, col = "lightgray")
  polygon(c(xval,rev(xval)),c(LL_all[[1]],rev(UL_all[[1]])),col=alpha(use_cols[1], 0.1),
          border = alpha(use_cols[1], 0.1))
  lines(x = xval, y = MID_all[[1]], col = use_cols[1])
  lines(Cases ~ Date, data = actual_dat_trim)
  legend("topleft", col = c(use_cols[1], "#000000"), lty = 1, 
         legend = c("Scenario 1", "Actual Cases"),
         bty = "n")
  dev.off()
}

make_graphs(forecast_fits, "forecast_fits.pdf")
make_graphs(forecast_fits_isolated70, "forecast_fits_isolated70.pdf")
  
### Estimated Deaths
# forecast
round( sapply( Map(comp_deaths, forecast_fits$UL), sum ), 0 )
round( sapply( Map(comp_deaths, forecast_fits$MID), sum ), 0 )
round( sapply( Map(comp_deaths, forecast_fits$LL), sum ), 0 )

# forecast with isolation
round( sapply( Map(comp_deaths, forecast_fits_isolated70$UL, forecaster = TRUE), sum ), 0 )
round( sapply( Map(comp_deaths, forecast_fits_isolated70$MID, forecaster = TRUE), sum ), 0 )
round( sapply( Map(comp_deaths, forecast_fits_isolated70$LL, forecaster = TRUE), sum ), 0 )

