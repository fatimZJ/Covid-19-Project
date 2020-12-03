## Load packages
library(deSolve)
library(DEoptim)
#library(tidyverse)
## Adapted from https://gist.github.com/jonocarroll/b17ce021b0637a31f584ed08a1fbe733
read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)[-1]
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}

jg_dat <- read.tcsv("data/dat_seir_code.csv")
jg_dat$Date <- as.Date(jg_dat$Date, format = "%a %d %b %Y") # Reformat date column

## Load the data
source('code/1_loadData.r')

acc_cases <- read.csv("data/CountData.csv")

## Load the get beta function

source("code/function_IEMAGSEIR_lsoda_scalars.R")


nlm_fun <- function(scalars, Acc_Cases) {
  
 # scalars <- exp(scalars)  
  
  test <- simulation_SEIR_model(R0t = 3.4,
                                dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                                dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                                dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                                dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                                POP = Irlpop,
                                numWeekStagger = c(3,6,9,12,15),
                                contacts_ireland = contacts,
                                dt = 1,  
                                tmax = 225,#length(Acc_Cases),
                                scalars = scalars)
  
  test_Cc <- rowSums(test$sol_out[grepl('Cc_',names(test$sol_out))])
  
  sum((Acc_Cases - (test_Cc[-1]))^2) 
  
}


scalars_init <- c(1.91365769, 1.16790530, 0.20215852, 
                  0.06962944, 0.18328977, 0.43466483,
                  0.60612880, 0.32053803, 0.45793943)

lower <- rep(0, length(scalars_init))
upper <- rep(3, length(scalars_init))
startt <- Sys.time()
ests_DEoptim <- DEoptim(nlm_fun, Acc_Cases = acc_cases$cases[1:225],
                        lower, upper, DEoptim.control(itermax = 200, reltol = 1e-8,
                                                      steptol = 10))

endt <- Sys.time()


nlm_fun(scalars = c(1.91365769, 1.16790530, 0.20215852, 
                    0.06962944, 0.18328977, 0.43466483,
                    0.60612880, 0.32053803, 0.45793943), Acc_Cases = acc_cases$cases[1:225])

plot(acc_cases$cases, type = "l")

ests_nlm <- nlm(nlm_fun,log(scalars_init),
                stepmax = 0.5, iterlim = 1000, Acc_Cases = acc_cases$cases[1:225])


## scalars nlm for 1:225 with different initial values
scalars <- c(0.322261050, 3.196447761, 0.210362425, 0.004546121, 0.240910420, 
             0.075327292,1.418424178, 0.267256734, 0.451187440)
# min 94888918

## scalars nlm for 1:225
scalars <- c(0.955845581, 2.144257251, 0.212053027, 
             0.064083360, 0.007385239, 0.292367212, 
             1.382529104, 0.269396475, 0.451376694)

## scalars nlm full range of acc_cases
scalars <- c(2.22010663, 0.94775863, 0.23377529, 
             0.07789591, 0.23403541, 0.43689101,
             0.47946681, 0.39151948, 0.36748479)

library(optimx)

startt <- Sys.time()
ests_optimx <-  optimx(c(1.000000, 0.97171108, 0.20812880, 0.07611468, 0.22169280, 
                         0.40394184, 0.44685362, 0.36467352, 0.42729580), 
                       Acc_Cases = acc_cases$cases[1:225],
                       nlm_fun, lower = 0, upper = Inf, control = list(maximize = FALSE))

endt <- Sys.time()

## DEoptim 1:225, upper = 3 reltol = 1e-6, steptol = 10
c(2.572162, 0.810708,  0.228309,   
0.021862, 0.505094, 0.134786,    
0.014523, 0.316367, 1.196971)
## DEoptim 1:225, upper = 3 reltol = 1e-4, steptol = 5

#Iteration: 21 bestvalit: 4243514669.109088 bestmemit:    1.574830    1.498398    0.159623    0.281354    0.214719    0.894793    0.015734    0.131745    0.380075
c(1.574830, 1.498398, 0.159623, 0.281354, 0.214719, 0.894793, 0.015734, 0.131745, 0.380075)

## DEoptim 1:225, upper = 1 reltol = 1e-4, steptol = 5

 scalars <- c(0.758034,  0.436740,  0.673103, 0.114594, 0.089738, 0.029895,
 0.384559, 0.650709, 0.467569)
 
test <- simulation_SEIR_model(R0t = 3.4,
                              dateStartSchoolClosure = as.Date('2020-03-12') , #schools closed before imtense lockdown
                              dateStartIntenseIntervention = as.Date('2020-03-27') , #Intense intervention: starts at Wuhan Lockdown
                              dateEndIntenseIntervention = as.Date('2020-05-18'), #date we begin relaxing intense intervention
                              dateStart = as.Date('2020-02-28'), #start date for epidemic in Ireland
                              POP = Irlpop,
                              numWeekStagger = c(3,6,9,12,15),
                              contacts_ireland = contacts,
                              dt = 1,  
                              tmax = 225,#length(acc_cases$cases),
                              scalars = c(2.572162, 0.810708,  0.228309,   
                                          0.021862, 0.505094, 0.134786,    
                                          0.014523, 0.316367, 1.196971)) 
                                #c(0.322261050, 3.196447761, 0.210362425, 0.004546121, 0.240910420, 
                                # 0.075327292,1.418424178, 0.267256734, 0.451187440))
                                #c(0.955845581, 2.144257251, 0.212053027, 
                                 # 0.064083360, 0.007385239, 0.292367212, 
                                  #1.382529104, 0.269396475, 0.451376694))#c(2.22010663, 0.94775863, 0.23377529, 
#0.07789591, 0.23403541, 0.43689101,
#0.47946681, 0.39151948, 0.36748479))

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
Cc <- test$sol_out[grepl('Cc_', names(test$sol_out))]

#plot(jg_dat$Susceptible[1:225], type = "l", lwd = 2)
#lines(test$sol_out$time, rowSums(S), type = "l", lwd = 2, col = "red")

plot(acc_cases$cases[1:225], type = "l", lwd = 2)
lines(test$sol_out$time, rowSums(Cc), type = "l", col = "red", lwd = 2)

plot(abs(acc_cases$cases[1:225] - rowSums(Cc)[-1]), type = "l", ylab = "cases", lwd = 2)
## Infected

plot(test$sol_out$time, (rowSums(cbind(Ip, IA, Ii, It, Iti, Iq))),lwd=2, type = "l", 
     xlab ="Time(days)",
     ylab = "Daily no. of infections", col='tomato',
     panel.first = rect(c(test$tStartSchoolClosure, test$tStartIntenseIntervention,
                          test$tEndIntenseIntervention, test$tRelaxIntervention1, 
                          test$tRelaxIntervention2,test$tRelaxIntervention3,
                          test$tRelaxIntervention4), -1e6,
                        c(test$tStartIntenseIntervention,test$tEndIntenseIntervention, 
                          test$tRelaxIntervention1,test$tRelaxIntervention2,
                          test$tRelaxIntervention3,test$tRelaxIntervention4,
                          test$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))

lines(jg_dat$Infected[1:225],lwd = 2)

## Susceptible

plot(test$sol_out$time, (rowSums(S)),lwd=2, type = "l", 
     xlab ="Time(days)",
     ylab = "Daily no. of susceptible", col='tomato',
     panel.first = rect(c(test$tStartSchoolClosure, test$tStartIntenseIntervention,
                          test$tEndIntenseIntervention, test$tRelaxIntervention1, 
                          test$tRelaxIntervention2,test$tRelaxIntervention3,
                          test$tRelaxIntervention4), -1e6,
                        c(test$tStartIntenseIntervention,test$tEndIntenseIntervention, 
                          test$tRelaxIntervention1,test$tRelaxIntervention2,
                          test$tRelaxIntervention3,test$tRelaxIntervention4,
                          test$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))

lines(jg_dat$Susceptible[1:225],lwd = 2)


## Recovered

plot(test$sol_out$time, (rowSums(R)),lwd=2, type = "l", 
     xlab ="Time(days)",
     ylab = "Daily no. of recovered", col='tomato',
     panel.first = rect(c(test$tStartSchoolClosure, test$tStartIntenseIntervention,
                          test$tEndIntenseIntervention, test$tRelaxIntervention1, 
                          test$tRelaxIntervention2,test$tRelaxIntervention3,
                          test$tRelaxIntervention4), -1e6,
                        c(test$tStartIntenseIntervention,test$tEndIntenseIntervention, 
                          test$tRelaxIntervention1,test$tRelaxIntervention2,
                          test$tRelaxIntervention3,test$tRelaxIntervention4,
                          test$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))

lines(jg_dat$Recovered[1:225],lwd = 2)

## Accumulated daily cases

plot(test$sol_out$time, (rowSums(Cc)),lwd=2, type = "l", 
     xlab ="Time(days)",
     ylab = "Daily cumulative no. of cases", col='tomato',
     panel.first = rect(c(test$tStartSchoolClosure, test$tStartIntenseIntervention,
                          test$tEndIntenseIntervention, test$tRelaxIntervention1, 
                          test$tRelaxIntervention2,test$tRelaxIntervention3,
                          test$tRelaxIntervention4), -1e6,
                        c(test$tStartIntenseIntervention,test$tEndIntenseIntervention, 
                          test$tRelaxIntervention1,test$tRelaxIntervention2,
                          test$tRelaxIntervention3,test$tRelaxIntervention4,
                          test$tRelaxIntervention5 ), 1e6, col=c('gray63','gray48', 
                                                                 'gray70', 'gray75',
                                                                 'gray80','gray85', 
                                                                 'gray90'), border=NA))
lines(cumsum(jg_dat$New.cases.confirmed[1:225]),lwd = 2)





