
## Load packages
library(deSolve)
library(plotly)
## Load the data
source('code/1_loadData.r')

## Load the get beta function
source("code/getbeta.R")

## Load the model function
source("code/function_SEIR_Model.R")

## load the simulation function 
source("code/function_SEIR_Simulation.R")

#scalars_test <- c(2.45199033, 0.920820963, 0.19089041, 0.01511209, 0.32126002, 
 #                 0.42580101, 0.52733763,0.30954179, 0.26381850)

## nlm results
#scalars_est <- c(2.27619341, 0.94894146, 0.21845530, 0.01503771, 0.31684429, 0.40807037, 0.49130460,
#  0.34598300, 0.20890798)

#scalars_est <- c(0.403832639,	2.644294363,	0.204475828,	0.001418374,	
#                 0.000557027,	0.724489601,	0.420480905,	0.35805272,	
#                 0.202294526)

#scalars_est <- c(2.4242360, 0.9130557, 0.2019896, 0.0153258, 0.3168725, 0.4115858,
#  0.5386053, 0.3075835, 0.2579594)
##NM Results 

#scalars_est <- c(2.59509019,0.996704393,	0.077085162,	0.751827556,
#                 0.018125642,	0.511875573,	0.450964092,	0.11619769,	0.065715592
#)
#scalars_est <- c(2.033784e+00, 1.093872e+00, 2.143144e-01, 3.064505e-02, 3.295509e-01, 3.654255e-01,
#  5.485281e-01, 3.371615e-01, 2.035889e-01)

## NM 5 seconds timeout
scalars_est <- c(1.652554584,	1.356798336, 0.202063943,	0.225574192,	0.004418746,	
  0.464799905,	0.552744885,	0.311931601,	0.258379996)

## NM 10 seconds timeout
scalars_est <- c(2.05442977,	1.054453831,	0.226260352,	0.003082808,	0.099852024,	0.574552131,	
  0.45266616,	0.349173341,	0.200236175)

names(scalars_est) <- unique(interventions_info$policy)


Base <- simulation_SEIR_model(R0t = 3.4,
                              POP = population,
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

#p

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

