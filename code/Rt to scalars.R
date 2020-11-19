## Load packages
library(deSolve)
#library(tidyverse)

## Load the data
source('code/1_loadData.r')

Beta <- read.csv("Data/JG_0.1_betas.csv")[,1]
JG_Rt <- read.csv("Data/JG_0.1_Rt.csv")[,2]

## Load the get beta function
source("code/getbeta.R")
source("Code/Test_fun.R")
source("Code/function_IEMAGSEIR_lsoda_scalars.R")

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

S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
Ev <- Base$sol_out[grepl('Ev_',names(Base$sol_out))]
Ip <- Base$sol_out[grepl('Ip_',names(Base$sol_out))]
IA <- Base$sol_out[grepl('IA_',names(Base$sol_out))]
Ii <- Base$sol_out[grepl('Ii_',names(Base$sol_out))]
It <- Base$sol_out[grepl('It_',names(Base$sol_out))]
Iti <- Base$sol_out[grepl('Iti_',names(Base$sol_out))]
Iq <- Base$sol_out[grepl('Iq_',names(Base$sol_out))]
R <- Base$sol_out[grepl('R_',names(Base$sol_out))]

plot(jg_dat$Infected[1:225], type = "l", lwd = 2, xlab ="Time(days)",
     ylab = "Daily no. of infections",
     panel.first = rect(c(Base$tStartSchoolClosure, Base$tStartIntenseIntervention,Base$tEndIntenseIntervention,
                          Base$tRelaxIntervention1, Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4), -1e6,
                        c(Base$tStartIntenseIntervention,Base$tEndIntenseIntervention, Base$tRelaxIntervention1,
                          Base$tRelaxIntervention2,Base$tRelaxIntervention3,
                          Base$tRelaxIntervention4,Base$tRelaxIntervention5 ), 1e6,
                        col=c('gray63','gray48', 'gray70', 'gray75','gray80','gray85', 'gray90'), border=NA))

lines(Base$sol_out$time[-1], (rowSums(cbind(Ip, IA, Ii,It,Iti,Iq)))[-1],lwd=2,col='tomato')
legend(120, 12000,legend = c("JG's model", "Our Model"),
       col = c("black", "tomato"), bty = 'n',lty = c(1,1),lwd = c(2,2), cex = 1)


################################################################################
nlm_fun <- function(scalars, Rt) {
  scalars <- exp(scalars)  
  ## Get base values given Rt 
  Base <- simulation_SEIR_model_ContRt(Rt = Rt,
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
  
  Base_S <- Base$sol_out[grepl('S_',names(Base$sol_out))]
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
 
 test_S <- test$sol_out[grepl('S_',names(test$sol_out))]
 
 norm((Base_S - test_S), type = "2")
  
}
#nlm_fun(scalars = log(c(1.10500000, 0.48585554, 0.1040644, 0.03805734,
#                0.11084640, 0.20197092, 0.22342681, 0.18233676,
#                0.21364790)), Rt = JG_Rt)

ests <- nlm(nlm_fun,log(c(1.10500000, 0.48585554, 0.1040644, 0.03805734,
                                        0.11084640, 0.20197092, 0.22342681, 0.18233676,
                                        0.21364790)),  stepmax = 0.5, Rt = JG_Rt)
            
    
