##### Estimate economic costs

setwd("ForecastApp/")
### Dates and costs from econ work
NI <- difftime( as.Date('2020-03-11'), as.Date('2020-02-29'), units = 'days' )[[1]]
L2 <- difftime( as.Date('2020-09-17'), as.Date('2020-08-18'), units = 'days' )[[1]]
L3 <- difftime( as.Date('2020-10-20'), as.Date('2020-09-18'), units = 'days' )[[1]]
L5 <- difftime( as.Date('2020-11-30'), as.Date('2020-10-21'), units = 'days' )[[1]]

days <- c(NI, L2, L3, L5) 
costs <- c(-2.47, 3.69, 0.9, 8.03)

daily_cost <- costs/days

L1_cost <- mean(daily_cost[1:2])
L4_cost <- mean(daily_cost[3:4])

policy <- c('No Intervention', paste0('Lockdown Level ', 1:5))
all_daily_costs <- c(daily_cost[1], L1_cost, daily_cost[2:3], 
                     L4_cost, daily_cost[4]) 

econ_df <- data.frame(Level = policy, Estimated_Costs = all_daily_costs)
write_csv(econ_df, 'data/Estimated_Costs.csv')

