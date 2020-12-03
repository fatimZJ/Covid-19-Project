## Load packages
library(deSolve)
#library(tidyverse)

## Load the data
source('code/1_loadData.r')

Beta <- read.csv("data/JG_0.1_betas.csv")[,1]
JG_Rt <- read.csv("data/JG_0.1_Rt.csv")[,1]

## Load the get beta function
source("code/getbeta.R")
source("code/Test_fun.R")
source("code/function_IEMAGSEIR_lsoda_scalars.R")

Base <- simulation_SEIR_model_ContRt(Rt = JG_Rt,
                                     dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                     dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                     dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                     dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                     POP = Irlpop,
                                     numWeekStagger = c(3,6,9,12,15),
                                     contacts_ireland = contacts,
                                     dt = 0.1,  
                                     tmax = 225,
                                     scalars = c(1, 1, 1, 1, 1, 1, 1, 1, 1)) 

## Plotting the solution

Base_S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
Base_Ev <- Base$sol_out[grepl('Ev_',names(Base$sol_out))]
Base_Ip <- Base$sol_out[grepl('Ip_',names(Base$sol_out))]
Base_IA <- Base$sol_out[grepl('IA_',names(Base$sol_out))]
Base_Ii <- Base$sol_out[grepl('Ii_',names(Base$sol_out))]
Base_It <- Base$sol_out[grepl('It_',names(Base$sol_out))]
Base_Iti <- Base$sol_out[grepl('Iti_',names(Base$sol_out))]
Base_Iq <- Base$sol_out[grepl('Iq_',names(Base$sol_out))]
Base_R <- Base$sol_out[grepl('R_',names(Base$sol_out))]

plot(jg_dat$Infected[1:225], type = "l", lwd = 2, xlab ="Time(days)",
     ylab = "Daily no. of infections",
     panel.first = rect(c(Base$tStartSchoolClosure, Base$tStartIntenseIntervention,Base$tEndIntenseIntervention,
                          Base$tRelaxIntervention1, Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4), -1e6,
                        c(Base$tStartIntenseIntervention,Base$tEndIntenseIntervention, Base$tRelaxIntervention1,
                          Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4,Base$tRelaxIntervention5 ), 1e6,
                        col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))

lines(Base$sol_out$time[-1], (rowSums(cbind(Base_Ip, Base_IA, Base_Ii, Base_It,
                                            Base_Iti,Base_Iq)))[-1],lwd=2,col='tomato')
legend(120, 12000,legend = c("JG's model", "Our Model"),
       col = c("black", "tomato"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)

base <- as.data.frame(cbind("time" = Base$sol_out$time, "Base_S" = rowSums(Base_S), 
                            "Base_Ev" = rowSums(Base_Ev), "Base_Ip" = rowSums(Base_Ip), 
                            "Base_IA" = rowSums(Base_IA), "Base_Ii" = rowSums(Base_Ii), 
                            "Base_It" = rowSums(Base_It), "Base_Iti" = rowSums(Base_Iti), 
                            "Base_Iq" = rowSums(Base_Iq), "Base_R" = rowSums(Base_R)))

plot_ly(data = base,  x = ~time, y = ~(rowSums(cbind(Base_Ip, Base_IA, Base_Ii, 
                                                     Base_It,Base_Iti,Base_Iq))),
        type = "scatter", mode = 'lines',  #height = 800,
        line = list(color = 'black'))%>%
  layout(yaxis = list(title = "Daily no. of infections"),
         xaxis = list(title ="Time(days)"))

################################################################################
## Get base values given Rt 
# Base <- simulation_SEIR_model_ContRt(Rt = Rt,
#                                      dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
#                                      dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
#                                      dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
#                                      dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
#                                      POP = Irlpop,
#                                      numWeekStagger = c(3,6,9,12,15),
#                                      contacts_ireland = contacts,
#                                      dt = 0.1,  
#                                      tmax = 225,
#                                      scalars = c(1, 1, 1, 1, 1, 1, 1, 1, 1)) 
# 
# Base_S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
# Base_Ev <- Base$sol_out[grepl('Ev_',names(Base$sol_out))]
# Base_Ip <- Base$sol_out[grepl('Ip_',names(Base$sol_out))]
# Base_IA <- Base$sol_out[grepl('IA_',names(Base$sol_out))]
# Base_Ii <- Base$sol_out[grepl('Ii_',names(Base$sol_out))]
# Base_It <- Base$sol_out[grepl('It_',names(Base$sol_out))]
# Base_Iti <- Base$sol_out[grepl('Iti_',names(Base$sol_out))]
# Base_Iq <- Base$sol_out[grepl('Iq_',names(Base$sol_out))]
# Base_R <- Base$sol_out[grepl('R_',names(Base$sol_out))]

nlm_fun <- function(scalars, Base) {
 # scalars <- exp(scalars)  
   
 # Base_S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
  ##
 test <- simulation_SEIR_model(R0t = 3.4,#JG_Rt[1],
                        dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                        dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                        dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                        dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                        POP = Irlpop,
                        numWeekStagger = c(3,6,9,12,15),
                        contacts_ireland = contacts,
                        dt = 0.1,  
                        tmax = 225,
                        scalars = scalars)
 
 test_S <- rowSums(test$sol_out[grepl('S_',names(test$sol_out))])
 
 norm((base$Base_S - test_S), type = "2")
  
}
#plot(test_S)
nlm_fun(scalars = c(1.91365769, 1.16790530, 0.20215852, 
                    0.06962944, 0.18328977, 0.43466483,
                    0.60612880, 0.32053803, 0.45793943), Base = Base)

#ests_nlm <- nlm(nlm_fun,log(c(2.21000000, 0.97171108, 0.20812880, 0.07611468, 0.22169280, 
#                          0.40394184, 0.44685362, 0.36467352, 0.42729580)),  
#            stepmax = 0.5, iterlim = 1000)

library(optimx)

ests_optimx <-  optimx(c(2.21000000, 0.97171108, 0.20812880, 0.07611468, 0.22169280, 
                         0.40394184, 0.44685362, 0.36467352, 0.42729580), Base = Base,
                       nlm_fun, lower = 0, upper = Inf, control = list(maximize = FALSE))
  
# p1        p2        p3 p4        p5        p6        p7
# L-BFGS-B 2.156007 0.9895324 0.2088416  0 0.3022948 0.4904968 0.5201892
# p8        p9   value fevals gevals niter convcode  kkt1
# L-BFGS-B 0.3302468 0.4579886 44607.1    178    178    NA        0 FALSE
# kkt2    xtime
# L-BFGS-B TRUE 15040.53
## nlm            
#c(1.91365769, 1.16790530, 0.20215852, 0.06962944, 0.18328977, 0.43466483,
# 0.60612880, 0.32053803, 0.45793943)

## Optimx

# p1       p2        p3         p4        p5        p6
#1.913658 1.167905 0.2021585 0.06962944 0.1832898 0.4346648 0.6061288 0.320538 0.4579394

## plotting the two 
plot(Base$sol_out$time, (rowSums(Base_S)),lwd=2, type = "l", xlab ="Time(days)",
     ylab = "Daily no. of infections",
     panel.first = rect(c(Base$tStartSchoolClosure, Base$tStartIntenseIntervention,
                          Base$tEndIntenseIntervention, Base$tRelaxIntervention1, 
                          Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4), -1e6,
                        c(Base$tStartIntenseIntervention,Base$tEndIntenseIntervention, 
                          Base$tRelaxIntervention1,Base$tRelaxIntervention2,
                          Base$tRelaxIntervention3,Base$tRelaxIntervention4,
                          Base$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))

lines(test$sol_out$time, (rowSums(S)),lwd=2,col='blue')


plot(Base$sol_out$time, (rowSums(cbind(Base_Ip, Base_IA, Base_Ii, Base_It,
                                       Base_Iti,Base_Iq))),lwd=2, type = "l", 
     xlab ="Time(days)",
     ylab = "Daily no. of infections",
     panel.first = rect(c(Base$tStartSchoolClosure, Base$tStartIntenseIntervention,
                          Base$tEndIntenseIntervention, Base$tRelaxIntervention1, 
                          Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4), -1e6,
                        c(Base$tStartIntenseIntervention,Base$tEndIntenseIntervention, 
                          Base$tRelaxIntervention1,Base$tRelaxIntervention2,
                          Base$tRelaxIntervention3,Base$tRelaxIntervention4,
                          Base$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))



test <- simulation_SEIR_model(R0t = 3.4,
                              dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                              dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                              dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                              dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                              POP = Irlpop,
                              numWeekStagger = c(3,6,9,12,15),
                              contacts_ireland = contacts,
                              dt = 0.1,  
                              tmax = 225,
                              scalars = c(1.91365769, 1.16790530, 0.20215852, 
                                          0.06962944, 0.18328977, 0.43466483,
                                          0.60612880, 0.32053803, 0.45793943))
                                

## Plotting the solution

S <- test$sol_out[grepl('S_',names(test$sol_out))]
Ev <- test$sol_out[grepl('Ev_',names(test$sol_out))]
Ip <- test$sol_out[grepl('Ip_',names(test$sol_out))]
IA <- test$sol_out[grepl('IA_',names(test$sol_out))]
Ii <- test$sol_out[grepl('Ii_',names(test$sol_out))]
It <- test$sol_out[grepl('It_',names(test$sol_out))]
Iti <- test$sol_out[grepl('Iti_',names(test$sol_out))]
Iq <- test$sol_out[grepl('Iq_',names(test$sol_out))]
R <- test$sol_out[grepl('R_',names(test$sol_out))]


plot(Base$sol_out$time, (rowSums(cbind(Base_Ip, Base_IA, Base_Ii, Base_It,
                                       Base_Iti,Base_Iq))),lwd=2, type = "l", 
     xlab ="Time(days)",
     ylab = "Daily no. of infections",
     panel.first = rect(c(Base$tStartSchoolClosure, Base$tStartIntenseIntervention,
                          Base$tEndIntenseIntervention, Base$tRelaxIntervention1, 
                          Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4), -1e6,
                        c(Base$tStartIntenseIntervention,Base$tEndIntenseIntervention, 
                          Base$tRelaxIntervention1,Base$tRelaxIntervention2,
                          Base$tRelaxIntervention3,Base$tRelaxIntervention4,
                          Base$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))

lines(test$sol_out$time, (rowSums(cbind(Ip, IA, Ii, It,
                                        Iti, Iq))),lwd=2,col='blue')

