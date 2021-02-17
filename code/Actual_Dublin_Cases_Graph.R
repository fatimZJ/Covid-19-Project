library(tidyverse)
dat <- read_csv("data/Daily_County_Cases.csv") %>% 
  filter(CountyName == "Dublin") %>%
  select(TimeStamp, ConfirmedCovidCases)
dat <- dat[-1, ]
dat$Date <- as.Date(substr(dat$TimeStamp, 1, 10), format = "%Y/%m/%d")
dat$Cases <- c(0, diff(dat$ConfirmedCovidCases))
dat_trim <- filter(dat, (Date >= as.Date("2020-12-01")) & (Date <= as.Date("2021-01-25")))

pdf("outputs/Actual_Dublin_Cases.pdf")
plot(dat_trim$Date, dat_trim$Cases, type = "l", main = "Actual Daily COVID-19 Cases Dublin",
     xlab = "Date", ylab = "Cases", ylim = c(0, max(dat_trim$Cases)),
     col = "blue")
grid()
dev.off()
