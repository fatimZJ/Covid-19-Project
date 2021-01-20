## Load packages
library(deSolve)
library(plotly)
library(Matrix)
library(tidyverse)

## Load the data
source('code/function_load_Data.r')
Data_list <- load_data()
list2env(Data_list, envir = .GlobalEnv)
rm(Data_list)
## Load the get beta function
source("code/function_getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
source("code/function_SEIR_Simulation.R")

## NM 10 seconds timeout

#scalars_est <- c(2.175648, 0.7734942, 0.2349296, 0.05953305, 0.1763019, 0.3924366, 0.4126519,
#                 0.4111283, 0.1557194)

scalars_est <- c(0.8645919, 1.757336, 0.2124244, 0.1005269, 0.1699474, 0.334918, 
                 0.5456985, 0.3506512, 0.2056358) 

names(scalars_est) <- unique(interventions_info$policy)

model_params <- c(3.6, 5.8, 13.0, 0.55, 0.05, 0.05, 0.21, 0.8, 0.1, 2)
names(model_params) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

Base <- simulation_SEIR_model(R0t = 3.4,
                              POP = population,
                              pars = model_params,
                              contacts_ireland = contacts,
                              interventions = interventions_info, 
                              dt = 1, 
                              dateStart = interventions_info$start[1],
                              dateEnd = interventions_info$end[dim(interventions_info)[1]],
                              scalars = scalars_est) 

## Plotting the solution

Cc <- Base$sol_out[grepl('Cc_',names(Base$sol_out))]

plot(cumulative_cases$cases, type = "l", lwd = 2, xlab ="Time(days)",
     ylab = "Daily no. of cumulative cases",
     panel.first = rect(c(Base$tstart_intervention[[2]], Base$tstart_intervention[[3]],
                          Base$tstart_intervention[[4]], Base$tstart_intervention[[5]], 
                          Base$tstart_intervention[[6]],Base$tstart_intervention[[7]],
                          Base$tstart_intervention[[8]], Base$tstart_intervention[[9]]), -1e6,
                        c(Base$tend_intervention[[2]]+1, Base$tend_intervention[[3]]+1,
                          Base$tend_intervention[[4]]+1, Base$tend_intervention[[5]]+1, 
                          Base$tend_intervention[[6]]+1,Base$tend_intervention[[7]]+1,
                          Base$tend_intervention[[8]]+1, Base$tend_intervention[[9]])+1, 1e6,
                        col=c('gray90','gray48', 'gray70', 'gray75','gray85','gray80', 'gray67','gray54'), border=NA))

lines(Base$sol_out$time[-1], (rowSums(Cc))[-1],lwd=2,col='tomato')
legend(40, 24000,legend = c("Actual data", "Our Model"),
       col = c("black", "red"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)

Sim <- as.data.frame(cbind("time" = Base$sol_out$time, 
                    "Cc" = rowSums(Cc)))

p <- plot_ly(data = Sim,  x = ~time, y = ~(Cc), name = 'Estimated',
             type = "scatter", mode = 'lines', # height = 480,
             line = list(color = 'tomato')) %>%
  layout(yaxis = list(title = "Daily cumulative no. of cases", titlefont = list(size = 18),
                      tickfont = list(size = 13)),
         xaxis = list(title ="Time(days)", titlefont = list(size = 18),
                      tickfont = list(size = 13))) %>% 
  add_trace(data = cumulative_cases, x = ~seq_along(date), y = ~cases, 
            name = 'Reported',
            mode = 'lines', 
            line = list(color = 'black') )

p <- p %>% layout(legend = list(x = 0.1, y = 0.9, font = list(size = 16)))

p

#plot(Base$sol_out$time, (rowSums(Cc)),lwd=2,col='tomato', type = "l")

layout(p, #title = 'Highlighting with Rectangles',
       shapes = list(
         list(type = "rect",
              fillcolor = "gray", line = list(color = "rgb(195, 195, 195)"), opacity = 0.15,
              x0 = Base$tstart_intervention[[2]], x1 = Base$tstart_intervention[[3]] -0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.5,
              x0 = Base$tstart_intervention[[3]]+0.1, x1 = Base$tstart_intervention[[4]] - 0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.4,
              x0 = Base$tstart_intervention[[4]] +0.1, x1 = Base$tstart_intervention[[5]]- 0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.3,
              x0 = Base$tstart_intervention[[5]]+ .1, x1 = Base$tstart_intervention[[6]]- 0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.25,
              x0 = Base$tstart_intervention[[6]]+ .1, x1 = Base$tstart_intervention[[7]]- 0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.35,
              x0 = Base$tstart_intervention[[7]]+ .1, x1 = Base$tstart_intervention[[8]]- 0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.44,
              x0 = Base$tstart_intervention[[8]]+ .1, x1 = Base$tstart_intervention[[9]]- 0.1, xref = "x",
              y0 = 0, y1 = 30000, yref = "y"),
         list(type = "rect",
              fillcolor = "gray", line = list(color = "gray"), opacity = 0.48,
              x0 = Base$tstart_intervention[[9]] + .1, x1 = Base$tend_intervention[[9]], xref = "x",
              y0 = 0, y1 = 30000, yref = "y")))#%>%
  #layout(xaxis = list(showgrid = F),
   #      yaxis = list(showgrid = F))

