
## Load the data
source('code/1_loadData.r')
contacts_ireland <- contacts
source("code/getbeta.R")
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

tmax <- 225


scalars <- c(0.955845581, 2.144257251, 0.212053027, 
  0.064083360, 0.007385239, 0.292367212, 
  1.382529104, 0.269396475, 0.451376694)
## Constraints

dateStartSchoolClosure = as.Date('2020-03-12') # School closure
dateStartIntenseIntervention = as.Date('2020-03-27')  #Intense intervention
dateEndIntenseIntervention = as.Date('2020-05-18') #date we begin relaxing intense intervention 
dateStart = as.Date('2020-02-28') #start date for epidemic in Ireland

numWeekStagger = c(3,6,9,12,15)

tStartSchoolClosure = as.vector(dateStartSchoolClosure - dateStart)+1 #Time point to add the school closure effect
tStartIntenseIntervention = as.vector(dateStartIntenseIntervention - dateStart)+1 # #Time point to add the intense lockdown effect
tEndIntenseIntervention = as.vector(dateEndIntenseIntervention - dateStart)+1     # for pw = 0.1
tRelaxIntervention1 = tEndIntenseIntervention + numWeekStagger[1]*7                               # for pw = 0.25
tRelaxIntervention2 = tEndIntenseIntervention + numWeekStagger[2]*7                               # for pw = 0.5
tRelaxIntervention3 = tEndIntenseIntervention + numWeekStagger[3]*7
tRelaxIntervention4 = tEndIntenseIntervention + numWeekStagger[4]*7                               # for pw = 0.25
tRelaxIntervention5 = tEndIntenseIntervention + numWeekStagger[5]*7                               # for pw = 0.5

## Model Parameters

pars <- c(4.9, 5.9, 7.0, 0.25, 0.05, 0.05, 0.5, 0.75, 0.13, 3.6)
names(pars) <- c("L","Cv","Dv","h","i","j","f","tv","q","TT")

Beta0 <- getbeta(R0t = 3.4, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)

Beta <- vector()
R <- vector()

## Get beta for each time point
for (t in 1:225){

  if(t < tStartSchoolClosure)
  {
    INTERVENTION <- scalars[1]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  # Schools closed before lockdown period
  if(t >= tStartSchoolClosure & t < tStartIntenseIntervention)
  {
    INTERVENTION <- scalars[2]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Intense intervention- lock down
  if(t >= tStartIntenseIntervention & t < tEndIntenseIntervention)
  {
    INTERVENTION <- scalars[3]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 1
  if(t >= tEndIntenseIntervention & t < tRelaxIntervention1)
  {
    INTERVENTION <- scalars[4]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 2
  if(t >= tRelaxIntervention1 & t < tRelaxIntervention2)
  {
    INTERVENTION <- scalars[5]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 3
  if(t >= tRelaxIntervention2 & t < tRelaxIntervention3)
  {
    INTERVENTION <- scalars[6]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 4
  if(t >= tRelaxIntervention3 & t < tRelaxIntervention4)
  {
    INTERVENTION <- scalars[7]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 5
  if(t >= tRelaxIntervention4 & t < tRelaxIntervention5)
  {
    INTERVENTION <- scalars[8]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Post lockdown, no intervention
  if(t >= tRelaxIntervention5)
  {
    INTERVENTION <- scalars[9]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  
  print(INTERVENTION)
  #print(t)
  C <- list()
  # total contacts matrix (work + school + household + other)
  C[[1]] <- CONSTRAINT%*%contacts_ireland[[1]]
  C[[2]] <- CONSTRAINT%*%contacts_ireland[[2]]
  C[[3]] <- CONSTRAINT%*%contacts_ireland[[3]]
  C[[4]] <- CONSTRAINT%*%contacts_ireland[[4]]
  C[[5]] <- CONSTRAINT%*%contacts_ireland[[5]]
  ## Estimating Beta
  R[t] <- getR0t(beta = Beta0, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = C)
  
  Beta[t] <- getbeta(R0t = 3.4, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = C)
}


plot(Beta, type = "l", lwd = 2)
plot(R, type = "l", lwd = 2)
#lines(jg_dat$Beta[1:tmax])

plot(jg_dat$Beta[1:tmax], type = "l", lwd = 2)
lines(Beta, col = "red", lwd = 2)

plot(jg_dat$Rt[1:tmax], type = "l", lwd = 2)
lines(R, col = "red", lwd = 2)
## Fitting a model with these beta values and no contacts matrices should give the 
## same output as fitting a model with a fixed beta value and these contact matrices

R0 <- 3.4

Beta1 <- vector()
R1 <- vector()

## Get beta for each time point
for (t in 1:225){
  
  if(t < tStartSchoolClosure)
  {
    INTERVENTION <- scalars[1]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  # Schools closed before lockdown period
  if(t >= tStartSchoolClosure & t < tStartIntenseIntervention)
  {
    INTERVENTION <- scalars[2]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Intense intervention- lock down
  if(t >= tStartIntenseIntervention & t < tEndIntenseIntervention)
  {
    INTERVENTION <- scalars[3]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 1
  if(t >= tEndIntenseIntervention & t < tRelaxIntervention1)
  {
    INTERVENTION <- scalars[4]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 2
  if(t >= tRelaxIntervention1 & t < tRelaxIntervention2)
  {
    INTERVENTION <- scalars[5]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 3
  if(t >= tRelaxIntervention2 & t < tRelaxIntervention3)
  {
    INTERVENTION <- scalars[6]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 4
  if(t >= tRelaxIntervention3 & t < tRelaxIntervention4)
  {
    INTERVENTION <- scalars[7]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Relaxing interventions - Phase 5
  if(t >= tRelaxIntervention4 & t < tRelaxIntervention5)
  {
    INTERVENTION <- scalars[8]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  #Post lockdown, no intervention
  if(t >= tRelaxIntervention5)
  {
    INTERVENTION <- scalars[9]
    CONSTRAINT  <- diag(INTERVENTION,16,16)
  }
  
  print(INTERVENTION)
  #print(t)
  
  ## Estimating Beta
  R1[t] <- R0 * INTERVENTION
  
  Beta1[t] <- Beta0 * INTERVENTION
}

Beta0 <- getbeta(R0t = R1, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = contacts_ireland)
plot(Beta1, type = "l", lwd = 2)
lines(Beta0, col = "tomato", lwd = 2)

R0 <- getR0t(beta = Beta0, pars = pars, p_age = Irlpop$propage, CONTACTMATRIX = C)
plot(R0, type = "l")
lines(R, col = "tomato")

plot(Beta, type = "l", lwd = 2)
lines(Beta1, col = "tomato", lwd = 2)

plot(R, type = "l", lwd = 2)
lines(R1, col = "tomato", lwd = 2)

