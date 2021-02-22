##### Create Graphs for Static Forecast Data for Manuscript

### Load in libraries and data
library(scales)
library(tidyverse)
actual_dat <- read_csv("data/Daily_County_Cases.csv") %>% 
  filter(CountyName == "Dublin") %>%
  select(TimeStamp, ConfirmedCovidCases)
setwd("Joint_Shiny_App/")
source("global.R")
load("data/forecast_fits.Rdata")
load("data/forecast_fits_isolated75.Rdata")
actual_dat <- actual_dat[-1, ]
actual_dat$Date <- as.Date(substr(actual_dat$TimeStamp, 1, 10), format = "%Y/%m/%d")
actual_dat$Cases <- c(0, diff(actual_dat$ConfirmedCovidCases))
actual_dat_trim <- filter(actual_dat, (Date >= as.Date("2020-12-01")) & (Date <= as.Date("2021-01-25")))

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
  
  ### Add in 'All' and 'Symptomatic' compartments
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
    x[["Sy"]][(nrow(x) - n_date):nrow(x)]
  }
  
  UL_y <- lapply(UL[-1], truncate_dat)
  MID_y <- lapply(MID[-1], truncate_dat)
  LL_y <- lapply(LL[-1], truncate_dat)
  
  max_ind <- which.max( sapply(UL_y, max) )
  
  pdf(filename)
  plot(x = xval, y = UL_y[[max_ind]], type = "n", ylab = "", xaxt = "n",
       xlab = "Date", ylim = c(0, max(UL_y[[max_ind]])), xlim = xval_disp,
       main = "Symptomatic Infections Forecast")
  axis.Date(1, at = xax, las = 2)
  abline(v = xax, lty = 2, col = "lightgray")
  abline(h = seq(0, 5000, 1000), lty = 2, col = "lightgray")
  
  # Add forecasts
  use_cols <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442")
  
  for (i in 1:4) {
    polygon(c(xval,rev(xval)),c(LL_y[[i]],rev(UL_y[[i]])),col=alpha(use_cols[i], 0.1),
          border = alpha(use_cols[i], 0.1))
  }
  
  for (i in 1:4) {
    lines(x = xval, y = MID_y[[i]], col = use_cols[i])
  }
  
  lines(Cases ~ Date, data = actual_dat_trim)
  
  legend("topleft", col = c(use_cols, "#000000"), lty = 1, 
         legend = c(paste("Scenario", 2:5), "Actual Cases"),
         bty = "n")
  dev.off()
}

make_graphs(forecast_fits, "forecast_fits.pdf")
make_graphs(forecast_fits_isolated75, "forecast_fits_isolated75.pdf")
  