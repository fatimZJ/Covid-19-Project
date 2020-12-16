
## Load packages
library(deSolve)

## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
source("code/function_SEIR_Simulation.R")



scalars_test <- c(2.45199033, 0.920820963, 0.19089041, 0.01511209, 0.32126002, 
                  0.42580101, 0.52733763,0.30954179, 0.26381850)


names(scalars_test) <- unique(interventions_info$policy)


Base <- simulation_SEIR_model(R0t = 3.4,
                              POP = population,
                              contacts_ireland = contacts,
                              interventions = interventions_info, 
                              dt = 1, 
                              dateStart = interventions_info$start[1],
                              dateEnd = interventions_info$end[dim(interventions_info)[1]],
                              scalars = scalars_test) 

## Plotting the solution

Cc <- Base$sol_out[grepl('Cc_',names(Base$sol_out))]

plot(cumulative_cases$cases, type = "l", lwd = 2, xlab ="Time(days)",
     ylab = "Daily no. of cumulative cases",
     panel.first = rect(c(Base$tstart_intervention[[2]], Base$tstart_intervention[[3]],
                          Base$tstart_intervention[[4]], Base$tstart_intervention[[5]], 
                          Base$tstart_intervention[[6]],Base$tstart_intervention[[7]],
                          Base$tstart_intervention[[8]], Base$tstart_intervention[[9]]), -1e6,
                        c(Base$tend_intervention[[2]], Base$tend_intervention[[3]],
                          Base$tend_intervention[[4]], Base$tend_intervention[[5]], 
                          Base$tend_intervention[[6]],Base$tend_intervention[[7]],
                          Base$tend_intervention[[8]], Base$tend_intervention[[9]]), 1e6,
                        col=c('gray88','gray48', 'gray70', 'gray75','gray80','gray85', 'gray70','gray54'), border=NA))

lines(Base$sol_out$time[-1], (rowSums(Cc))[-1],lwd=2,col='tomato')
legend(40, 24000,legend = c("Actual data", "Our Model"),
       col = c("black", "red"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)



#plot(Base$sol_out$time, (rowSums(Cc)),lwd=2,col='tomato', type = "l")



