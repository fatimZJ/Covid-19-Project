## Reproducing paper plots
par(mfrow = c(1,1))
load('outputs/SEIcIscR/AGEcovid_IDurInf7.rdata')


length(AGEcovid_IDurInf7sc)

names(AGEcovid_IDurInf7sc[[1]]$summary)

doNothing <- (AGEcovid_IDurInf7sc[[1]]$summary$median)

base <- (AGEcovid_IDurInf7sc[[2]]$summary$median)

march <- (AGEcovid_IDurInf7sc[[3]]$summary$median)

april <- (AGEcovid_IDurInf7sc[[4]]$summary$median)

# incidence over time
agegp =1

plot(doNothing[,agegp], type='l', lwd=2,
     main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
     xlab="Time(days)", ylab="Daily no. of infections");
lines(base[,agegp],lwd=2,col='grey40')
lines(march[,agegp],lwd=2,col='steelblue')
lines(april[,agegp],lwd=2,col='tomato',lty='dashed')


load('outputs/SEIcIscR/covid_DurInf7.rdata')


length(covid_DurInf7sc)

dim(covid_DurInf7sc[[1]]$inc$med)


doNothing <- (covid_DurInf7sc[[1]]$inc$med)

base <- (covid_DurInf7sc[[2]]$inc$med)

march <- (covid_DurInf7sc[[3]]$inc$med)

april <- (covid_DurInf7sc[[4]]$inc$med)




load('outputs/SEIcIscR/covid_DurInf7I2000.rdata')


length(covid_DurInf7I2000sc)

dim(covid_DurInf7I2000sc[[1]]$inc$med)


doNothing <- (covid_DurInf7I2000sc[[1]]$inc$med)

base <- (covid_DurInf7I2000sc[[2]]$inc$med)

march <- (covid_DurInf7I2000sc[[3]]$inc$med)

april <- (covid_DurInf7I2000sc[[4]]$inc$med)


# incidence over time
agegp = 13

plot(doNothing[,agegp], type='l', lwd=2,
     main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
     xlab="Time(days)", ylab="Daily no. of infections");
lines(base[,agegp],lwd=2,col='grey40')
lines(march[,agegp],lwd=2,col='steelblue')
lines(april[,agegp],lwd=2,col='tomato',lty='dashed')


load('outputs/SEIcIscR/covid_IDurInf7.rdata')

length(covid_IDurInf7sc)

covid_IDurInf7sc[[1]]$summary$median
plot(covid_IDurInf7sc[[1]]$summary$median/100, type = "l")


plot(ovid_IDurInf7sc[[1]]$summary$median/100, type='l', lwd=2,
     xlab="Time(days)", ylab="Daily no. of infections");
lines(covid_IDurInf7sc[[2]]$summary$median/100,lwd=2,col='grey40')
lines(covid_IDurInf7sc[[3]]$summary$median/100,lwd=2,col='steelblue')
lines(covid_IDurInf7sc[[4]]$summary$median/100,lwd=2,col='tomato',lty='dashed')

load('outputs/SEIcIscR/covid_IDurInf7I2000.rdata')

length(covid_IDurInf7I2000sc)

covid_IDurInf7I2000sc[[1]]$summary$median
plot(covid_IDurInf7I2000sc[[1]]$summary$median/100, type = "l")


plot(covid_IDurInf7I2000sc[[1]]$summary$median/100, type='l', lwd=2,
     xlab="Time(days)", ylab="Daily no. of infections");
lines(covid_IDurInf7I2000sc[[2]]$summary$median/100,lwd=2,col='grey40')
lines(covid_IDurInf7I2000sc[[3]]$summary$median/100,lwd=2,col='steelblue')
lines(covid_IDurInf7I2000sc[[4]]$summary$median/100,lwd=2,col='tomato',lty='dashed')

################################################################################

load('outputs/SEIR/covid_DurInf7.rdata')

par(mfrow=c(1,1))
length(covid_DurInf7)

dim(covid_DurInf7[[1]]$inc$med)


doNothing <- (covid_DurInf7[[1]]$inc$med)

base <- (covid_DurInf7[[2]]$inc$med)

march <- (covid_DurInf7[[3]]$inc$med)

april <- (covid_DurInf7[[4]]$inc$med)

# incidence over time
agegp = 1
plot(doNothing[,agegp], type='l', lwd=2,
     main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
     xlab="Time(days)", ylab="Daily no. of infections");
lines(base[,agegp],lwd=2,col='grey40')
lines(march[,agegp],lwd=2,col='steelblue')
lines(april[,agegp],lwd=2,col='tomato',lty='dashed')


doNothing <- (covid_DurInf7[[1]]$inc$med)

base <- (covid_DurInf7[[2]]$inc$med)

march <- (covid_DurInf7[[3]]$inc$med)

april <- (covid_DurInf7[[4]]$inc$med)

# incidence over time
agegp = c(4,5)


plot(rowSums(doNothing[,agegp]/1000), type='l', lwd=2,
     main=paste0("Incidence for age ", (agegp[1]-1)*5," to <", agegp[length(agegp)]*5),
     xlab="Time(days)", ylab="New cases per day (in thousands)", ylim = c(0,1));
#     panel.first = rect(c(1,7), -1e6, c(3,10), 1e6, col='green', border=NA));
lines(rowSums(base[,agegp]/1000),lwd=2,col='grey40')
lines(rowSums(march[,agegp]/1000),lwd=2,col='steelblue')
lines(rowSums(april[,agegp]/1000),lwd=2,col='tomato',lty='dashed')



load('outputs/SEIcIscR/covid_DurInf7.rdata')


length(covid_DurInf7sc)

dim(covid_DurInf7sc[[1]]$inc$med)


doNothing <- (covid_DurInf7sc[[1]]$inc$med)

base <- (covid_DurInf7sc[[2]]$inc$med)

march <- (covid_DurInf7sc[[3]]$inc$med)

april <- (covid_DurInf7sc[[4]]$inc$med)




load('outputs/SEIR/covid_DurInf7I2000.rdata')


length(covid_DurInf7I2000)

dim(covid_DurInf7I2000[[1]]$inc$med)


doNothing <- (covid_DurInf7I2000[[1]]$inc$med)

base <- (covid_DurInf7I2000[[2]]$inc$med)

march <- (covid_DurInf7I2000[[3]]$inc$med)

april <- (covid_DurInf7I2000[[4]]$inc$med)


# incidence over time
agegp = 6
plot(doNothing[,agegp], type='l', lwd=2,
     main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
     xlab="Time(days)", ylab="Daily no. of infections");
lines(base[,agegp],lwd=2,col='grey40')
lines(march[,agegp],lwd=2,col='steelblue')
lines(april[,agegp],lwd=2,col='tomato',lty='dashed')

#### Infection period 7 days
load('outputs/SEIR/covid_DurInf7.rdata')

par(mfrow=c(1,1))
length(covid_DurInf7)

dim(covid_DurInf7[[1]]$inc$med)

doNothing <- (covid_DurInf7[[1]]$inc$med)

base <- (covid_DurInf7[[2]]$inc$med)

march <- (covid_DurInf7[[3]]$inc$med)

april <- (covid_DurInf7[[4]]$inc$med)

# incidence over time
agegp = 1
plot(doNothing[,agegp], type='l', lwd=2,
     main=paste0("Incidence for age [",(agegp-1)*5,',',agegp*5,')'),
     xlab="Time(days)", ylab="Daily no. of infections");
lines(base[,agegp],lwd=2,col='grey40')
lines(march[,agegp],lwd=2,col='steelblue')
lines(april[,agegp],lwd=2,col='tomato',lty='dashed')

# grouped
agegp = c(4,5)

plot(rowSums(doNothing[,agegp]/1000), type='l', lwd=2,
     main=paste0("Incidence for age ", (agegp[1]-1)*5," to <", agegp[length(agegp)]*5),
     xlab="Time(days)", ylab="New cases per day (in thousands)", ylim = c(0,1));
#     panel.first = rect(c(1,7), -1e6, c(3,10), 1e6, col='green', border=NA));
lines(rowSums(base[,agegp]/1000),lwd=2,col='grey40')
lines(rowSums(march[,agegp]/1000),lwd=2,col='steelblue')
lines(rowSums(april[,agegp]/1000),lwd=2,col='tomato',lty='dashed')

